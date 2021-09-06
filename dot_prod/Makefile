###################################################################
#  Makefile for Example: dot-product
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


cpp : dot_prod_cpp.exe dot_prod_mpi_simple_cpp.exe \
      dot_prod_mpi_cpp.exe dot_prod_mpi_fancy_cpp.exe \
      dot_prod_mpi_noreduce_cpp.exe

dot_prod_cpp.exe : dot_prod.cpp get_time.c
	$(CXX) $(CXXFLAGS) dot_prod.cpp get_time.c -lm -o $@

dot_prod_mpi_simple_cpp.exe : dot_prod_mpi_simple.cpp
	$(MPICXX) $(CXXFLAGS) dot_prod_mpi_simple.cpp -lm -o $@

dot_prod_mpi_cpp.exe : dot_prod_mpi.cpp
	$(MPICXX) $(CXXFLAGS) dot_prod_mpi.cpp -lm -o $@

dot_prod_mpi_fancy_cpp.exe : dot_prod_mpi_fancy.cpp
	$(MPICXX) $(CXXFLAGS) dot_prod_mpi_fancy.cpp -lm -o $@

dot_prod_mpi_noreduce_cpp.exe : dot_prod_mpi_noreduce.cpp
	$(MPICXX) $(CXXFLAGS) dot_prod_mpi_noreduce.cpp -lm -o $@


c : dot_prod_c.exe dot_prod_mpi_simple_c.exe \
    dot_prod_mpi_c.exe dot_prod_mpi_fancy_c.exe \
    dot_prod_mpi_noreduce_c.exe

dot_prod_c.exe : dot_prod.c get_time.c
	$(CC) $(CFLAGS) dot_prod.c get_time.c -lm -o $@

dot_prod_mpi_simple_c.exe : dot_prod_mpi_simple.c
	$(MPICC) $(CFLAGS) dot_prod_mpi_simple.c -lm -o $@

dot_prod_mpi_c.exe : dot_prod_mpi.c
	$(MPICC) $(CFLAGS) dot_prod_mpi.c -lm -o $@

dot_prod_mpi_fancy_c.exe : dot_prod_mpi_fancy.c
	$(MPICC) $(CFLAGS) dot_prod_mpi_fancy.c -lm -o $@

dot_prod_mpi_noreduce_c.exe : dot_prod_mpi_noreduce.c
	$(MPICC) $(CFLAGS) dot_prod_mpi_noreduce.c -lm -o $@


f90 : dot_prod_f90.exe dot_prod_mpi_simple_f90.exe \
      dot_prod_mpi_f90.exe dot_prod_mpi_fancy_f90.exe 

dot_prod_f90.exe : dot_prod.f90 get_time.f90
	$(F90) $(FFLAGS) dot_prod.f90 get_time.f90 -o $@

dot_prod_mpi_simple_f90.exe : dot_prod_mpi_simple.f90
	$(MPIF90) $(FFLAGS) dot_prod_mpi_simple.f90 -o $@

dot_prod_mpi_f90.exe : dot_prod_mpi.f90
	$(MPIF90) $(FFLAGS) dot_prod_mpi.f90 -o $@

dot_prod_mpi_fancy_f90.exe : dot_prod_mpi_fancy.f90
	$(MPIF90) $(FFLAGS) dot_prod_mpi_fancy.f90 -o $@


clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
