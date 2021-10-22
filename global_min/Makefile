###################################################################
#  Makefile for Example: global minimization
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 4370/6370
#  7 February 2015
###################################################################

# compilers and flags
CXX = g++
CXXFLAGS = -O2


# makefile targets
all : glob_min_cpp.exe


glob_min_cpp.exe : glob_min.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@


clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
