###################################################################
#  Makefile for Example: computing pi
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 6370
#  26 February 2013
###################################################################

F90 = gfortran
MPIF90 = mpif90
FFLAGS = -O2

CXX = g++
MPICXX = mpicxx
CXXFLAGS = -O2

CC = gcc
MPICC = mpicc
CFLAGS = -O2


# makefile targets
all : cpp c f90


cpp : pi_comp_cpp.exe pi_comp_mpi_simple_cpp.exe pi_comp_mpi_cpp.exe

pi_comp_cpp.exe : pi_comp.cpp get_time.c
	$(CXX) $(CXXFLAGS) pi_comp.cpp get_time.c -lm -o $@

pi_comp_mpi_simple_cpp.exe : pi_comp_mpi_simple.cpp
	$(MPICXX) $(CXXFLAGS) pi_comp_mpi_simple.cpp -o $@

pi_comp_mpi_cpp.exe : pi_comp_mpi.cpp
	$(MPICXX) $(CXXFLAGS) pi_comp_mpi.cpp -o $@


c : pi_comp_c.exe pi_comp_mpi_simple_c.exe pi_comp_mpi_c.exe

pi_comp_c.exe : pi_comp.c get_time.c
	$(CC) $(CFLAGS) pi_comp.c get_time.c -lm -o $@

pi_comp_mpi_simple_c.exe : pi_comp_mpi_simple.c
	$(MPICC) $(CFLAGS) pi_comp_mpi_simple.c -o $@

pi_comp_mpi_c.exe : pi_comp_mpi.c
	$(MPICC) $(CFLAGS) pi_comp_mpi.c -o $@


f90 : pi_comp_f90.exe pi_comp_mpi_simple_f90.exe pi_comp_mpi_f90.exe

pi_comp_f90.exe : pi_comp.f90 get_time.f90
	$(F90) $(FFLAGS) pi_comp.f90 get_time.f90 -o $@

pi_comp_mpi_simple_f90.exe : pi_comp_mpi_simple.f90
	$(MPIF90) $(FFLAGS) pi_comp_mpi_simple.f90 -o $@

pi_comp_mpi_f90.exe : pi_comp_mpi.f90
	$(MPIF90) $(FFLAGS) pi_comp_mpi.f90 -o $@


clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
