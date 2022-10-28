#!/bin/bash
#
# File to set the relevant environment variables to compile/run this example on Maneframe2

module load intel/18.0.2 impi/18.0.2 hypre/2.19
export HYPRE_ROOT=$TACC_HYPRE_DIR/skylake
export MPICXX=mpicxx
