#!/bin/bash

# set desired vector lengths
N=(20000 200000 2000000 20000000)

# configure build
cmake .

# build executables
make

# loop over vector lengths
for n in "${N[@]}"
do

   # run serial test
   echo "  "
   echo "running serial version, dot_prod.serial, N = $n:"
   ./dot_prod.serial $n

   # run CUDA test
   echo "  "
   echo "running CUDA version, dot_prod.cuda, N = $n:"
   ./dot_prod.cuda $n
done