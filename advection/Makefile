###################################################################
#  Makefile for Example: 2D wave equation
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 4370/6370
#  7 February 2015
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


cpp : advection_cpp.exe advection_mpi_cpp.exe advection_1sided.exe

advection_cpp.exe : advection.cpp output.cpp initialize.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

advection_mpi_cpp.exe : advection_mpi.cpp output_mpi.cpp initialize_mpi.cpp check_err.c
	$(MPICXX) $(CXXFLAGS) $^ -lm -o $@

advection_1sided.exe : advection_1sided.cpp output_1sided.cpp initialize_1sided.cpp exchange_1sided.cpp check_err.c
	$(MPICXX) $(CXXFLAGS) $^ -lm -o $@


c : advection_c.exe advection_mpi_c.exe

advection_c.exe : advection.c output.c initialize.c get_time.c
	$(CC) $(CFLAGS) $^ -lm -o $@

advection_mpi_c.exe : advection_mpi.c output_mpi.c initialize_mpi.c check_err.c
	$(MPICC) $(CFLAGS) $^ -lm -o $@


f90 : advection_f90.exe advection_mpi_f90.exe

advection_f90.exe : advection.f90 output.f90 initialize.f90 get_time.f90
	$(F90) $(FFLAGS) $^ -o $@

advection_mpi_f90.exe : exchange.f90 advection_mpi.f90 initialize_mpi.f90 output_mpi.f90
	$(MPIF90) $(FFLAGS) $^ -o $@


clean :
	\rm -f *.o u_sol* *.mod *.png
	\rm -rf results_*

realclean : clean
	\rm -f *.exe *~ 


####### End of Makefile #######
