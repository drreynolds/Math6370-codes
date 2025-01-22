#!/bin/bash

# set desired number of chemical bins
N=10000

# configure build
cmake .

# build executables
make

# run serial test
echo "  "
echo "running serial version, chemistry.serial:"
./chemistry.serial $N

# run CUDA test
echo "  "
echo "running CUDA version, chemistry.cuda:"
./chemistry.cuda $N

# run UVM test
echo "  "
echo "running UVM version, chemistry.cuda:"
./chemistry.uvm $N
