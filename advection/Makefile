###################################################################
#  Makefile for Example: 2D wave equation
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 4370/6370
#  7 February 2015
###################################################################

# set compilers based on machine
ifeq ($(HOST), asymptote)
  CXX = g++-11
  CC  = gcc-11
  F90 = gfortran-11
else
  CXX = g++
  CC  = gcc
  F90 = gfortran
endif

# flags
CXXFLAGS = -O2
CFLAGS   = -O2
FFLAGS   = -O2

OMPFLAGS = -fopenmp


# makefile targets
all : advection_omp_cpp.exe advection_omp_c.exe advection_omp_f90.exe

advection_omp_cpp.exe : advection_omp.cpp output.cpp initialize.cpp get_time.c
	$(CXX) $(OMPFLAGS) $(CXXFLAGS) $^ -lm -o $@

advection_omp_c.exe : advection_omp.c output.c initialize.c get_time.c
	$(CC) $(OMPFLAGS) $(CFLAGS) $^ -lm -o $@

advection_omp_f90.exe : advection_omp.f90 output.f90 initialize.f90 get_time.f90
	$(F90) $(OMPFLAGS) $(FFLAGS) $^ -o $@


clean :
	\rm -f *.o u_sol*

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
