###################################################################
#  Makefile for Example: matrix-vector product
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 6370
#  26 February 2015
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

all_cpp : matvec_row_cpp.exe matvec_col_cpp.exe \
      matvec_row_auto_cpp.exe matvec_col_auto_cpp.exe \
      matvec_omp_cpp.exe

all_c : matvec_row_c.exe matvec_col_c.exe \
      matvec_row_auto_c.exe matvec_col_auto_c.exe \
      matvec_omp_c.exe

all_f90 : matvec_row_f90.exe matvec_col_f90.exe \
      matvec_row_auto_f90.exe matvec_col_auto_f90.exe \
      matvec_omp_f90.exe


matvec_row_cpp.exe : matvec_row.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

matvec_col_cpp.exe : matvec_col.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

matvec_row_auto_cpp.exe : matvec_row.cpp get_time.c
	$(CXX) $(AUTOPAR) $(CXXFLAGS) $^ -lm -o $@

matvec_col_auto_cpp.exe : matvec_col.cpp get_time.c
	$(CXX) $(AUTOPAR) $(CXXFLAGS) $^ -lm -o $@

matvec_omp_cpp.exe : matvec_omp.cpp get_time.c
	$(CXX) $(OMPFLAGS) $(CXXFLAGS) $^ -lm -o $@


matvec_row_c.exe : matvec_row.c get_time.c
	$(CC) $(CFLAGS) $^ -lm -o $@

matvec_col_c.exe : matvec_col.c get_time.c
	$(CC) $(CFLAGS) $^ -lm -o $@

matvec_row_auto_c.exe : matvec_row.c get_time.c
	$(CC) $(AUTOPAR) $(CFLAGS) $^ -lm -o $@

matvec_col_auto_c.exe : matvec_col.c get_time.c
	$(CC) $(AUTOPAR) $(CFLAGS) $^ -lm -o $@

matvec_omp_c.exe : matvec_omp.c get_time.c
	$(CC) $(OMPFLAGS) $(CFLAGS) $^ -lm -o $@


matvec_row_f90.exe : matvec_row.f90 get_time.f90
	$(F90) $(FFLAGS) $^ -o $@

matvec_col_f90.exe : matvec_col.f90 get_time.f90
	$(F90) $(FFLAGS) $^ -o $@

matvec_row_auto_f90.exe : matvec_row.f90 get_time.f90
	$(F90) $(AUTOPAR) $(FFLAGS) $^ -o $@

matvec_col_auto_f90.exe : matvec_col.f90 get_time.f90
	$(F90) $(AUTOPAR) $(FFLAGS) $^ -o $@

matvec_omp_f90.exe : matvec_omp.f90
	$(F90) $(OMPFLAGS) $(FFLAGS) $^ -o $@


clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~ b*.txt


####### End of Makefile #######
