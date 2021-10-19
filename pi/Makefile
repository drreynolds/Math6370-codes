###################################################################
#  Makefile for Example: computing pi
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 6370
#  4 January 2013
###################################################################

# compilers
CXX = g++

# flags
CXXFLAGS = -O2

# makefile targets
all : pi_comp_cpp.exe

pi_comp_cpp.exe : pi_comp.cpp get_time.c
	$(CXX) $(CXXFLAGS) pi_comp.cpp get_time.c -lm -o $@

clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
