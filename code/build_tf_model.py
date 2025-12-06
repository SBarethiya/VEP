""" builds the tensorflow execution graph for the model by parsing network_specs """
import os, warnings
warnings.simplefilter(action='ignore', category=DeprecationWarning)
from os.path import isfile
import argparse
import yaml
import tensorflow as tf
import numpy as np
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

import gen_structure_graph as gsg

def initializer(init, shape):
    if init == "zero":
        return tf.zeros(shape)

    elif init == "he":
        fan_in = np.prod(shape)
        std = 1 / np.sqrt(fan_in)
        return tf.random_uniform(shape, minval=-std, maxval=std, dtype=tf.float32)


def node_average_gc(inputs, adj_mtx, activation, filters=None, trainable=True):
    """
    Implementation of graph convolution operation, modified slightly from the implementation used in
    Protein Interface Prediction using Graph Convolutional Networks
    https://github.com/fouticus/pipgcn
    available under the MIT License, Copyright 2020 Alex Fout
    param:
        inputs: input tensor of shape (batch_size, number_of_vertices, encoding_len)
        adj_mtx: adjacency matrix tensor of shape (number_of_vertices, number_of_vertices)
        activation: activation function to apply after convolution
        filters: number of output filters
        trainable: whether the weights are trainable
    return:
        output_signal: output tensor of shape (batch_size, number_of_vertices, filters)
    """
    # node_average_gc_dist_thresh

    vertices = inputs  # shape: (batch_size, number_of_vertices, encoding_len)
    v_shape = vertices.get_shape()

    # create new weights # (v_dims, filters)
    center_weights = tf.Variable(initializer("he", (v_shape[-1].value, filters)), name="Wc", trainable=trainable)
    neighbor_weights = tf.Variable(initializer("he", (v_shape[-1].value, filters)), name="Wn", trainable=trainable)
    bias = tf.Variable(initializer("zero", (filters,)), name="b", trainable=trainable)

    # center signals are simply the center node value times the weight
    # shape: (batch_size, number_of_vertices, num_filters)
    center_signals = tf.reshape(tf.matmul(tf.reshape(vertices, (-1, v_shape[-1])),
                                          center_weights),
                                (-1, v_shape[1], filters))

    # apply neighbor weight to each neighbor
    # shape: (batch_size, number_of_vertices, num_filters)
    neighbor_signals_sep = tf.reshape(tf.matmul(tf.reshape(vertices, (-1, v_shape[-1])), neighbor_weights),
                                      (-1, v_shape[1], filters))

    # compute full neighbor signals
    neighbor_signals = tf.divide(tf.matmul(tf.tile(adj_mtx[None], (tf.shape(vertices)[0], 1, 1)),
                                           neighbor_signals_sep),
                                 tf.reshape(tf.maximum(tf.constant(1, dtype=tf.float32),
                                                       tf.reduce_sum(adj_mtx, axis=1)), (-1, 1)))

    # final output signal
    output_signal = activation(center_signals + neighbor_signals + bias)

    return output_signal

def eval_spec(spec, local_scope):
    """
    Recursively evaluates a specification dictionary, replacing any string that starts with "~" with the result of
    evaluating the string as a Python expression in the given local scope.
    param: 
        spec: specification dictionary, list, or string
        local_scope: dictionary representing the local scope for evaluation
    return: 
        evaluated specification
    """
    if isinstance(spec, dict):
        return {k: eval_spec(v, local_scope) for k, v in spec.items()}
    elif isinstance(spec, list):
        return [eval_spec(i, local_scope) for i in spec]
    elif isinstance(spec, str) and spec.startswith("~"):
        return eval(spec[1:], globals(), local_scope)
    return spec


def bg_inference(net_fn, adj_mtx, ph_inputs_dict):
    """
    builds the inference part of the graph by parsing the network specification yaml file
    param:
        net_fn: filename of the network specification yaml file
        args: dictionary of arguments (not used here, but could be useful for custom layer functions)
        adj_mtx: adjacency matrix tensor for graph convolutional layers (can be None if no graph layers are used)
        ph_inputs_dict: dictionary of placeholder input tensors
    return:
        predictions: output tensor of the network
    """
    with open(net_fn, "r") as f:
        yml = yaml.safe_load(f)

    layers = [ph_inputs_dict["raw_seqs"]]
    for layer_spec in yml["network"]:
        if layer_spec["layer_func"] == "~node_average_gc" and adj_mtx is None:
            raise ValueError("must specify a protein structure graph (adj_mtx) when using a graph convolutional layer")
        
        parsed_spec = eval_spec(layer_spec, locals())
        layer = parsed_spec["layer_func"](**parsed_spec["arguments"])
        layers.append(layer)
        
    layers.append(tf.layers.dense(layers[-1], units=1, activation=None, name="output"))
    predictions = tf.squeeze(layers[-1], axis=1)
    return predictions


def bg_loss(inf_graph):
    """ builds the graph by adding the required loss ops """
    scores = inf_graph["ph_inputs_dict"]["scores"]
    predictions = inf_graph["predictions"]
    loss = tf.compat.v1.losses.mean_squared_error(labels=scores,
                                        predictions=predictions,
                                        reduction=tf.compat.v1.losses.Reduction.MEAN)
    return loss


def bg_training(args, loss):
    """ adds operations needed for training """
    learning_rate = args["learning_rate"]

    # create the adam descent optimizer with the given learning rate
    optimizer = tf.compat.v1.train.AdamOptimizer(learning_rate)

    # create a variable to track the global step
    global_step = tf.Variable(0, name='global_step', trainable=False)

    # Use the optimizer to apply the gradients that minimize the loss
    train_op = optimizer.minimize(loss, global_step=global_step)

    return global_step, train_op


def get_placeholder_inputs(args,data_shape):
    """ sets up the place holder input variables, into which the actual data will be fed """
    # remove the batch dimension from the data shape (this can be inferred real time and is not needed for
    # constructing the net
    raw_seqs_shape = tuple([None] + list(data_shape[1:]))

    # put all placeholders in a dictionary
    ph_dict = {
        # placeholder for main input data          
        "raw_seqs": tf.compat.v1.placeholder(tf.float32, shape=raw_seqs_shape, name="raw_seqs_placeholder"),
        # the score of the example (enrichment ratio, score from Enrich2)
        "scores": tf.compat.v1.placeholder(tf.float32, shape=None, name="scores_placeholder"),
        # placeholder for dropout
        "training": tf.compat.v1.placeholder_with_default(False, shape=(), name="training_ph")            
    }

#    ros_score_shape = tuple([None] + list(ros_data_shape[1:])) 
#    ph_dict.update({"ros_score": tf.compat.v1.placeholder(tf.float32, shape=ros_score_shape, name="ros_score_placeholder")})
        
    return ph_dict


def get_placeholder_metrics():
    """ add placeholder for evaluation metrics to graph. these are used as a simple way to get the evaluation metrics
        added to the graph and tensorboard """
    # metrics computed outside of tensorflow and added to the graph for the tensorboard visualization
    mse_ph = tf.compat.v1.placeholder(tf.float32, name="mean_squared_error")
    pearsonr_ph = tf.compat.v1.placeholder(tf.float32, name="pearsonr")
    r2_ph = tf.compat.v1.placeholder(tf.float32, name="r2")
    spearmanr_ph = tf.compat.v1.placeholder(tf.float32, name="spearmanr")

    metrics_ph_dict = {"mse": mse_ph,
                       "pearsonr": pearsonr_ph,
                       "r2": r2_ph,
                       "spearmanr": spearmanr_ph}

    validation_loss_ph = tf.compat.v1.placeholder(tf.float32, name="validation_loss_placeholder")
    training_loss_ph = tf.compat.v1.placeholder(tf.float32, name="training_loss_placeholder")

    return metrics_ph_dict, validation_loss_ph, training_loss_ph


def bg_summaries(metrics_ph_dict, validation_loss_ph, training_loss_ph):
    """ add summary scalars to the graph for keeping track of various stats for tensorboard """
    # validation and training loss, evaluated at every epoch
    tf.compat.v1.summary.scalar("validation_loss", validation_loss_ph, collections=["summaries_per_epoch"])
    tf.compat.v1.summary.scalar("training_loss", training_loss_ph, collections=["summaries_per_epoch"])

    # metrics (not necessarily evaluated at every epoch, although I have been doing that)
    tf.compat.v1.summary.scalar("mse", metrics_ph_dict["mse"], collections=["summaries_metrics"])
    tf.compat.v1.summary.scalar("pearsonr", metrics_ph_dict["pearsonr"], collections=["summaries_metrics"])
    tf.compat.v1.summary.scalar("r2", metrics_ph_dict["r2"], collections=["summaries_metrics"])
    tf.compat.v1.summary.scalar("spearmanr", metrics_ph_dict["spearmanr"], collections=["summaries_metrics"])

    # build the summary Tensor based on the TF collection of summaries
    # this is used as the "op" to feed in values for the above metrics
    summaries_per_epoch = tf.compat.v1.summary.merge_all("summaries_per_epoch")
    summaries_metrics = tf.compat.v1.summary.merge_all("summaries_metrics")

    return summaries_per_epoch, summaries_metrics


def build_inference_graph(args, encoded_data_shape):
    """ builds the inference part of the graph. the encoded data shape is expected to have the first dimension
        be the number of examples (batch size). it will be ignored, but still expected, so make it 1 if needed """

    # load adjacency matrix for gcn
    graph_fn = args["graph_fn"]
    adj_mtx = None
    if isfile(graph_fn):
        g = gsg.load_graph(graph_fn)
        adj_mtx = gsg.ordered_adjacency_matrix(g)
        adj_mtx = tf.convert_to_tensor(adj_mtx)

    # placeholder for inputs (raw sequences, labels, weights, is_training, etc)
    ph_inputs_dict = get_placeholder_inputs(args, encoded_data_shape)
    
    # build the inference part of the graph that gets the output values from the inputs
    predictions = bg_inference(args["net_file"], adj_mtx, ph_inputs_dict)

    # place all relevant tensorflow variables and ops into a dictionary for passing around to other functions
    inf_graph = {"ph_inputs_dict": ph_inputs_dict,
                 "predictions": predictions}

    return inf_graph


def build_training_graph(args, inf_graph):
    """ builds the training part of the graph """

    # placeholders for metrics and validation, training loss
    metrics_ph_dict, validation_loss_ph, training_loss_ph = get_placeholder_metrics()
    summaries_per_epoch, summaries_metrics = bg_summaries(metrics_ph_dict, validation_loss_ph, training_loss_ph)

    # op for loss calculation
    loss = bg_loss(inf_graph)

    # ops that calculate and apply gradients as well as global step
    global_step, train_op = bg_training(args, loss)

    # variable initializer op for initializing network weights
    init_global = tf.compat.v1.global_variables_initializer()

    train_graph = {"loss": loss,
                   "global_step": global_step,
                   "train_op": train_op,
                   "init_global": init_global,
                   "summaries_per_epoch": summaries_per_epoch,
                   "summaries_metrics": summaries_metrics,
                   "validation_loss_ph": validation_loss_ph,
                   "training_loss_ph": training_loss_ph,
                   "metrics_ph_dict": metrics_ph_dict}

    return train_graph


def build_graph_from_args_dict(args, encoded_data_shape, reset_graph=True):
    """
    builds the tensorflow execution graph for the model by parsing network_specs
    param:
        args: dictionary of arguments
        encoded_data_shape: shape of the encoded data (including batch dimension)
        reset_graph: whether to reset the default graph before building
    return:
        inf_graph: inference part of the graph
        train_graph: training part of the graph
    """
    if isinstance(args, argparse.Namespace):
        args = vars(args)

    if reset_graph:
        tf.compat.v1.reset_default_graph()

    inf_graph = build_inference_graph(args, encoded_data_shape)
    train_graph = build_training_graph(args, inf_graph)
    return inf_graph, train_graph


def main():
    pass


if __name__ == "__main__":
    main()
