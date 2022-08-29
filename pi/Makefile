###################################################################
#  Makefile for Example: computing pi
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 4370 / 6370
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
all : pi_comp_cpp.exe pi_comp_c.exe pi_comp_f90.exe

pi_comp_cpp.exe : pi_comp.cpp get_time.c
	$(CXX) $(CXXFLAGS) pi_comp.cpp get_time.c -lm -o $@

pi_comp_c.exe : pi_comp.c get_time.c
	$(CC) $(CFLAGS) pi_comp.c get_time.c -lm -o $@

pi_comp_f90.exe : pi_comp.f90 get_time.f90
	$(F90) $(FFLAGS) pi_comp.f90 get_time.f90 -o $@

clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~

####### End of Makefile #######
