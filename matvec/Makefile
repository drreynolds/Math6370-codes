###################################################################
#  Makefile for Example: matrix-vector product
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 4370 / 6370
###################################################################

# F90
F90 = gfortran
MPIF90 = mpif90
FFLAGS = -O2

# C++
CXX = g++
MPICXX = mpicxx
CXXFLAGS = -O2

# C
CC = gcc
MPICC = mpicc
CFLAGS = -O2


# makefile targets
all : cpp c f90


cpp : matvec_row_cpp.exe matvec_col_cpp.exe \
      matvec_mpi1_cpp.exe matvec_mpi2_cpp.exe \
      matvec_mpi3_cpp.exe

matvec_row_cpp.exe : matvec_row.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

matvec_col_cpp.exe : matvec_col.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

matvec_mpi1_cpp.exe : matvec_mpi1.cpp
	$(MPICXX) $(CXXFLAGS) $^ -lm -o $@

matvec_mpi2_cpp.exe : matvec_mpi2.cpp
	$(MPICXX) $(CXXFLAGS) $^ -lm -o $@

matvec_mpi3_cpp.exe : matvec_mpi3.cpp
	$(MPICXX) $(CXXFLAGS) $^ -lm -o $@


c : matvec_row_c.exe matvec_col_c.exe \
    matvec_mpi1_c.exe matvec_mpi2_c.exe \
    matvec_mpi3_c.exe

matvec_row_c.exe : matvec_row.c get_time.c
	$(CC) $(CFLAGS) $^ -lm -o $@

matvec_col_c.exe : matvec_col.c get_time.c
	$(CC) $(CFLAGS) $^ -lm -o $@

matvec_mpi1_c.exe : matvec_mpi1.c
	$(MPICC) $(CFLAGS) $^ -lm -o $@

matvec_mpi2_c.exe : matvec_mpi2.c
	$(MPICC) $(CFLAGS) $^ -lm -o $@

matvec_mpi3_c.exe : matvec_mpi3.c
	$(MPICC) $(CFLAGS) $^ -lm -o $@


f90 : matvec_row_f90.exe matvec_col_f90.exe \
      matvec_mpi1_f90.exe matvec_mpi2_f90.exe \
      matvec_mpi3_f90.exe

matvec_col_f90.exe : matvec_col.f90 get_time.f90
	$(F90) $(FFLAGS) $^ -o $@

matvec_row_f90.exe : matvec_row.f90 get_time.f90
	$(F90) $(FFLAGS) $^ -o $@

matvec_mpi1_f90.exe : matvec_mpi1.f90
	$(MPIF90) $(FFLAGS) $^ -o $@

matvec_mpi2_f90.exe : matvec_mpi2.f90
	$(MPIF90) $(FFLAGS) $^ -o $@

matvec_mpi3_f90.exe : matvec_mpi3.f90
	$(MPIF90) $(FFLAGS) $^ -o $@


clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~ b*.txt


####### End of Makefile #######
