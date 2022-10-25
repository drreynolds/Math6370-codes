###################################################################
#  Makefile for Example: simple MPI program
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 6370
#  5 February 2013
###################################################################

F90 = gfortran
MPIF90 = mpif90
FFLAGS = -O2

CC = gcc
MPICC = mpicc
CFLAGS = -O2

CXX = g++
MPICXX = mpicxx
CXXFLAGS = -O2


# makefile targets
all : simple_cpp.exe simple_c.exe simple_f90.exe

simple_cpp.exe : simple.cpp
	$(MPICXX) $(CXXFLAGS) $^ -lm -o $@

simple_c.exe : simple.c
	$(MPICC) $(CFLAGS) $^ -lm -o $@

simple_f90.exe : simple.f90
	$(MPIF90) $(FFLAGS) $^ -o $@

clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
