###################################################################
#  Makefile for Example: computing chemical equlibria
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
CXXFLAGS = -O2 -std=c++11
CFLAGS   = -O2
FFLAGS   = -O2

# makefile targets
all : chemistry_cpp.exe chemistry_lambda.exe chemistry_c.exe chemistry_f90.exe

chemistry_cpp.exe : chemistry.cpp chem_solver.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

chemistry_lambda.exe : chemistry.cpp chem_solver_lambda.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

chemistry_c.exe : chemistry.c chem_solver.c get_time.c
	$(CC) $(CFLAGS) $^ -lm -o $@

chemistry_f90.exe : chemistry.f90 chem_solver.f90 get_time.f90
	$(F90) $(FFLAGS) $^ -o $@

clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
