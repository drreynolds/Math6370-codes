#!/bin/bash

# set desired vector lengths
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
   echo "running serial version, axpy.serial, N = $n:"
   ./axpy.serial $n

   # run CUDA test
   echo "  "
   echo "running CUDA version, axpy.cuda, N = $n:"
   ./axpy.cuda $n
done