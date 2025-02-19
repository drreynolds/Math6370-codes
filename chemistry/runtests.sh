#!/bin/bash

# set desired numbers of chemical bins
N=(10 100 1000 10000)

# configure build
cmake .

# build executables
make

# loop over vector lengths
for n in "${N[@]}"
do

   # run serial test
   echo "  "
   echo "running serial version, chemistry.serial, N = $n:"
   ./chemistry.serial $n

   # run CUDA test
   echo "  "
   echo "running CUDA version, chemistry.cuda, N = $n:"
   ./chemistry.cuda $n

   # run UVM test
   echo "  "
   echo "running UVM version, chemistry.uvm, N = $n:"
   ./chemistry.uvm $n
done