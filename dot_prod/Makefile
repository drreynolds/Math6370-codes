###################################################################
#  Makefile for Example: dot-product
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


all_cpp : dot_prod_cpp.exe dot_prod_auto_cpp.exe \
      dot_prod_omp_cpp.exe dot_prod_omp_simple_cpp.exe \
      dot_prod_omp2_cpp.exe

all_c : dot_prod_c.exe dot_prod_auto_c.exe \
      dot_prod_omp_c.exe dot_prod_omp_simple_c.exe \
      dot_prod_omp2_c.exe

all_f90 : dot_prod_f90.exe dot_prod_auto_f90.exe \
      dot_prod_omp_f90.exe dot_prod_omp_simple_f90.exe \
      dot_prod_omp2_f90.exe


dot_prod_cpp.exe : dot_prod.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

dot_prod_auto_cpp.exe : dot_prod.cpp get_time.c
	$(CXX) $(AUTOPAR) $(CXXFLAGS) $^ -lm -o $@

dot_prod_omp_cpp.exe : dot_prod_omp.cpp get_time.c
	$(CXX) $(OMPFLAGS) $(CXXFLAGS) $^ -lm -o $@

dot_prod_omp_simple_cpp.exe : dot_prod_omp_simple.cpp get_time.c
	$(CXX) $(OMPFLAGS) $(CXXFLAGS) $^ -lm -o $@

dot_prod_omp2_cpp.exe : dot_prod_omp2.cpp get_time.c
	$(CXX) $(OMPFLAGS) $(CXXFLAGS) $^ -lm -o $@


dot_prod_c.exe : dot_prod.c get_time.c
	$(CC) $(CFLAGS) $^ -lm -o $@

dot_prod_auto_c.exe : dot_prod.c get_time.c
	$(CC) $(AUTOPAR) $(CFLAGS) $^ -lm -o $@

dot_prod_omp_c.exe : dot_prod_omp.c get_time.c
	$(CC) $(OMPFLAGS) $(CFLAGS) $^ -lm -o $@

dot_prod_omp_simple_c.exe : dot_prod_omp_simple.c get_time.c
	$(CC) $(OMPFLAGS) $(CFLAGS) $^ -lm -o $@

dot_prod_omp2_c.exe : dot_prod_omp2.c get_time.c
	$(CC) $(OMPFLAGS) $(CFLAGS) $^ -lm -o $@


dot_prod_f90.exe : dot_prod.f90 get_time.f90
	$(F90) $(FFLAGS) $^ -o $@

dot_prod_auto_f90.exe : dot_prod.f90 get_time.f90
	$(F90) $(AUTOPAR) $(FFLAGS) $^ -o $@

dot_prod_omp_f90.exe : dot_prod_omp.f90 get_time.f90
	$(F90) $(OMPFLAGS) $(FFLAGS) $^ -o $@

dot_prod_omp_simple_f90.exe : dot_prod_omp_simple.f90 get_time.f90
	$(F90) $(OMPFLAGS) $(FFLAGS) $^ -o $@

dot_prod_omp2_f90.exe : dot_prod_omp2.f90 get_time.f90
	$(F90) $(OMPFLAGS) $(FFLAGS) $^ -o $@


clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
