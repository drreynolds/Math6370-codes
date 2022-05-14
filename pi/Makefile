###################################################################
#  Makefile for Example: computing pi
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
AUTOPAR  = -floop-parallelize-all -ftree-parallelize-loops=4


# makefile targets
all : all_cpp all_c all_f90

all_cpp : pi_comp_cpp.exe pi_comp_auto_cpp.exe \
      pi_comp_omp_simple_cpp.exe pi_comp_omp_cpp.exe \
      pi_comp_omp_crit_cpp.exe

all_c : pi_comp_c.exe pi_comp_auto_c.exe \
      pi_comp_omp_simple_c.exe pi_comp_omp_c.exe \
      pi_comp_omp_crit_c.exe

all_f90 : pi_comp_f90.exe pi_comp_auto_f90.exe \
      pi_comp_omp_simple_f90.exe pi_comp_omp_f90.exe \
      pi_comp_omp_crit_f90.exe


pi_comp_cpp.exe : pi_comp.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

pi_comp_auto_cpp.exe : pi_comp.cpp get_time.c
	$(CXX) $(AUTOPAR) $(CXXFLAGS) $^ -lm -o $@

pi_comp_omp_simple_cpp.exe : pi_comp_omp_simple.cpp get_time.c
	$(CXX) $(OMPFLAGS) $(CXXFLAGS) $^ -lm -o $@

pi_comp_omp_cpp.exe : pi_comp_omp.cpp get_time.c
	$(CXX) $(OMPFLAGS) $(CXXFLAGS) $^ -lm -o $@

pi_comp_omp_crit_cpp.exe : pi_comp_omp_crit.cpp get_time.c
	$(CXX) $(OMPFLAGS) $(CXXFLAGS) $^ -lm -o $@


pi_comp_c.exe : pi_comp.c get_time.c
	$(CC) $(CFLAGS) $^ -lm -o $@

pi_comp_auto_c.exe : pi_comp.c get_time.c
	$(CC) $(AUTOPAR) $(CFLAGS) $^ -lm -o $@

pi_comp_omp_simple_c.exe : pi_comp_omp_simple.c get_time.c
	$(CC) $(OMPFLAGS) $(CFLAGS) $^ -lm -o $@

pi_comp_omp_c.exe : pi_comp_omp.c get_time.c
	$(CC) $(OMPFLAGS) $(CFLAGS) $^ -lm -o $@

pi_comp_omp_crit_c.exe : pi_comp_omp_crit.c get_time.c
	$(CC) $(OMPFLAGS) $(CFLAGS) $^ -lm -o $@


pi_comp_f90.exe : pi_comp.f90 get_time.f90
	$(F90) $(FFLAGS) $^ -o $@

pi_comp_auto_f90.exe : pi_comp.f90 get_time.f90
	$(F90) $(AUTOPAR) $(FFLAGS) $^ -o $@

pi_comp_omp_simple_f90.exe : pi_comp_omp_simple.f90 get_time.f90
	$(F90) $(OMPFLAGS) $(FFLAGS) $^ -o $@

pi_comp_omp_f90.exe : pi_comp_omp.f90 get_time.f90
	$(F90) $(OMPFLAGS) $(FFLAGS) $^ -o $@

pi_comp_omp_crit_f90.exe : pi_comp_omp_crit.f90 get_time.f90
	$(F90) $(OMPFLAGS) $(FFLAGS) $^ -o $@


clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
