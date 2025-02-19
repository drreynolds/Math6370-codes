#!/bin/bash

# set desired numbers of intervals
N=(10000 100000 1000000 10000000 100000000)

# configure build
cmake .

# build executables
make

# loop over vector lengths
for n in "${N[@]}"
do

   # run serial test
   echo "  "
   echo "running serial version, pi_comp.serial, N = $n:"
   ./pi_comp.serial $n

   # run CUDA test
   echo "  "
   echo "running CUDA version, pi_comp.cuda, N = $n:"
   ./pi_comp.cuda $n

done