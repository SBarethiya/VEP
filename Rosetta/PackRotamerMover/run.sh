#!/bin/sh
# Provide the last residue number
end=248
# Provide the file name
# filename="relaxed_200_0.pdb"
filename="../../data/avgfp/avgfp_rosetta_model.pdb"
# You can change this value if want to use more cpus
increment=10
# Change this value if you do not want to start from residue 1
start=1
int_end=$increment

while [ $start -le $end ]; do
  python ddG_calculation.py $start $int_end $filename
  echo $start $int_end $filename
  start=$((start + increment))
  int_end=$((int_end + increment))
done
