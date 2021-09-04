###################################################################
#  Makefile for Example: vector operations
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 4370/6370
#  14 January 2015
###################################################################

# compilers
CXX = g++
CC  = gcc
F90 = gfortran

# flags
CXXFLAGS = -O2
CFLAGS   = -O2
FFLAGS   = -O2

OMPFLAGS = -fopenmp
AUTOPAR  = -floop-parallelize-all -ftree-parallelize-loops=4


# makefile targets
all : all_cpp all_c all_f90

all_cpp : axpy_cpp.exe axpy_auto_cpp.exe axpy_omp_cpp.exe

all_c : axpy_c.exe axpy_auto_c.exe axpy_omp_c.exe

all_f90 : axpy_f90.exe axpy_auto_f90.exe axpy_omp_f90.exe


axpy_cpp.exe : axpy.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

axpy_auto_cpp.exe : axpy.cpp get_time.c
	$(CXX) $(CXXFLAGS) $(AUTOPAR) $^ -lm -o $@

axpy_omp_cpp.exe : axpy_omp.cpp get_time.c
	$(CXX) $(OMPFLAGS) $(CXXFLAGS) $^ -lm -o $@


axpy_c.exe : axpy.c get_time.c
	$(CC) $(CFLAGS) $^ -lm -o $@

axpy_auto_c.exe : axpy.c get_time.c
	$(CC) $(CFLAGS) $(AUTOPAR) $^ -lm -o $@

axpy_omp_c.exe : axpy_omp.c get_time.c
	$(CC) $(OMPFLAGS) $(CFLAGS) $^ -lm -o $@


axpy_f90.exe : axpy.f90 get_time.f90
	$(F90) $(FFLAGS) $^ -o $@

axpy_auto_f90.exe : axpy.f90 get_time.f90
	$(F90) $(FFLAGS) $(AUTOPAR) $^ -o $@

axpy_omp_f90.exe : axpy_omp.f90 get_time.f90
	$(F90) $(OMPFLAGS) $(FFLAGS) $^ -o $@


clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
