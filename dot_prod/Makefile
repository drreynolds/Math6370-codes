###################################################################
#  Makefile for Example: dot-product
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 6370
#  4 January 2013
###################################################################

# compilers
CXX = g++
CC  = gcc
F90 = gfortran

# flags 
CXXFLAGS = -O2
CFLAGS   = -O2
FFLAGS   = -O2

# makefile targets
all : dot_prod_cpp.exe dot_prod_c.exe dot_prod_f90.exe dot_prod_f.exe

dot_prod_cpp.exe : dot_prod.cpp get_time.c
	$(CXX) $(CXXFLAGS) dot_prod.cpp get_time.c -o $@

dot_prod_c.exe : dot_prod.c get_time.c
	$(CC) $(CFLAGS) dot_prod.c get_time.c -lm -o $@

dot_prod_f90.exe : dot_prod.f90 get_time.f90
	$(F90) $(FFLAGS) dot_prod.f90 get_time.f90 -o $@

dot_prod_f.exe : dot_prod.f get_time.f90
	$(F90) $(FFLAGS) dot_prod.f get_time.f90 -o $@

clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
