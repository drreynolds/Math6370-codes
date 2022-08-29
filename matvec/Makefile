###################################################################
#  Makefile for Example: matrix-vector product
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 4370 / 6370
###################################################################

# F90
F90 = gfortran
FFLAGS = -O2

# C++
CXX = g++
CXXFLAGS = -O2

# C
CC = gcc
CFLAGS = -O2


# makefile targets
all : cpp c f90

cpp : matvec_row_cpp.exe matvec_col_cpp.exe

matvec_row_cpp.exe : matvec_row.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

matvec_col_cpp.exe : matvec_col.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@



c : matvec_row_c.exe matvec_col_c.exe

matvec_row_c.exe : matvec_row.c get_time.c
	$(CC) $(CFLAGS) $^ -lm -o $@

matvec_col_c.exe : matvec_col.c get_time.c
	$(CC) $(CFLAGS) $^ -lm -o $@



f90 : matvec_row_f90.exe matvec_col_f90.exe

matvec_col_f90.exe : matvec_col.f90 get_time.f90
	$(F90) $(FFLAGS) $^ -o $@

matvec_row_f90.exe : matvec_row.f90 get_time.f90
	$(F90) $(FFLAGS) $^ -o $@



clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~ b*.txt


####### End of Makefile #######
