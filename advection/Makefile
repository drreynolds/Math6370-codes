###################################################################
#  Makefile for Example: 2D wave equation
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 4370/6370
#  18 October 2021
###################################################################

# compilers/flags
CXX = g++
CXXFLAGS = -O2
ifndef RAJA
  RAJA = /usr/local/raja-0.12.1
endif

# makefile targets
all : advection_cpp.exe

advection_cpp.exe : advection.cpp output.cpp initialize.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

clean :
	\rm -f *.o u_sol*

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
