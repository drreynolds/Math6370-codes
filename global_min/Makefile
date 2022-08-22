###################################################################
#  Makefile for Example: global minimization
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

OMPFLAGS = -fopenmp
AUTOPAR  = -floop-parallelize-all -ftree-parallelize-loops=4


# makefile targets
all : all_cpp all_c all_f90

all_cpp : glob_min_cpp.exe glob_min_auto_cpp.exe \
      glob_min2_cpp.exe glob_min_omp_cpp.exe

all_c : glob_min_c.exe glob_min_auto_c.exe \
      glob_min2_c.exe glob_min_omp_c.exe

all_f90 : glob_min_f90.exe glob_min_auto_f90.exe \
      glob_min2_f90.exe glob_min_omp_f90.exe


glob_min_cpp.exe : glob_min.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

glob_min_auto_cpp.exe : glob_min.cpp get_time.c
	$(CXX) $(CXXFLAGS) $(AUTOPAR) $^ -lm -o $@

glob_min2_cpp.exe : glob_min2.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

glob_min_omp_cpp.exe : glob_min_omp.cpp get_time.c
	$(CXX) $(OMPFLAGS) $(CXXFLAGS) $^ -lm -o $@


glob_min_c.exe : glob_min.c get_time.c
	$(CC) $(CFLAGS) $^ -lm -o $@

glob_min_auto_c.exe : glob_min.c get_time.c
	$(CC) $(CFLAGS) $(AUTOPAR) $^ -lm -o $@

glob_min2_c.exe : glob_min2.c get_time.c
	$(CC) $(CFLAGS) $^ -lm -o $@

glob_min_omp_c.exe : glob_min_omp.c get_time.c
	$(CC) $(OMPFLAGS) $(CFLAGS) $^ -lm -o $@


glob_min_f90.exe : glob_min.f90 get_time.f90
	$(F90) $(FFLAGS) $^ -o $@

glob_min_auto_f90.exe : glob_min.f90 get_time.f90
	$(F90) $(FFLAGS) $(AUTOPAR) $^ -o $@

glob_min2_f90.exe : glob_min2.f90 get_time.f90
	$(F90) $(FFLAGS) $^ -o $@

glob_min_omp_f90.exe : glob_min_omp.f90 get_time.f90
	$(F90) $(OMPFLAGS) $(FFLAGS) $^ -o $@


clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
