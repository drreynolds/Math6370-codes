###################################################################
#  Makefile for Example: 2D wave equation
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 6370
#  4 January 2013
###################################################################

# compilers/flags
CXX = g++
CC  = gcc
F90 = gfortran

# compilers/flags
CXXFLAGS = -O2
CFLAGS   = -O2
FFLAGS   = -O2


# makefile targets
all : advection_cpp.exe advection_c.exe advection_f90.exe

advection_cpp.exe : advection.cpp output.cpp initialize.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

advection_c.exe : advection.c output.c initialize.c get_time.c
	$(CC) $(CFLAGS) $^ -lm -o $@

advection_f90.exe : advection.f90 output.f90 initialize.f90 get_time.f90
	$(F90) $(FFLAGS) $^ -o $@

clean :
	\rm -f *.o u_sol*

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
