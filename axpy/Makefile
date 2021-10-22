###################################################################
#  Makefile for Example: vector operations
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 4370/6370
#  14 January 2015
###################################################################

# compilers and flags
CXX = g++
CXXFLAGS = -O2


# makefile targets
all : axpy_cpp.exe


axpy_cpp.exe : axpy.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@



clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
