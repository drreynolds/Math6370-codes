#!/bin/bash
#
# File to set the relevant environment variables to compile/run this example on Maneframe2

module load gcc-9.2 openmpi/3.1.6-liesbt2 hypre/2.20.0-4vrndhi
export HYPRE_ROOT=.
export MPICXX=mpicxx
