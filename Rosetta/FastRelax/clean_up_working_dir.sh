#!/bin/bash


cd ../rosetta_working_dir
rm -rf !(templates) 
rm -f structure_*.pdb
rm -f *mutate.xml
rm -f *mutation.resfile
rm -f *energy.txt
rm -f *score.sc
