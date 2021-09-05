###################################################################
#  Makefile for Example: openmp schedule example
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 4370/6370
#  7 February 2015
###################################################################

# C++
CXX = g++
CXXFLAGS = -O2

# C
CC = gcc
CFLAGS = -O2

# F90
F90 = gfortran
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
