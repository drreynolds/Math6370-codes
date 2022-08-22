###################################################################
#  Makefile for Example: openmp schedule example
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 4370/6370
#  7 February 2015
###################################################################

# set compilers based on machine
ifeq ($(HOST), asymptote)
  CXX = g++-12
  CC  = gcc-12
  F90 = gfortran-12
else
  CXX = g++
  CC  = gcc
  F90 = gfortran
endif

# flags
CXXFLAGS = -O2
CFLAGS   = -O2
FFLAGS   = -O2

# parallelization
OMPFLAGS = -fopenmp


# makefile targets
all : omp_schedule_cxx.exe omp_schedule_c.exe omp_schedule_f90.exe

omp_schedule_cxx.exe : omp_schedule.cpp
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $^ -lm -o $@

omp_schedule_c.exe : omp_schedule.c
	$(CC) $(CFLAGS) $(OMPFLAGS) $^ -lm -o $@

omp_schedule_f90.exe : omp_schedule.f90
	$(F90) $(FFLAGS) $(OMPFLAGS) $^ -lm -o $@


clean :
	\rm -f *.o

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
