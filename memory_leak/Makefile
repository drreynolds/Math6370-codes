###################################################################
#  Makefile for serial example: memory leak
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 4370/6370
#  23 January 2015
###################################################################

# C++
CXX = g++
CXXFLAGS = -O0 -p -g


# makefile targets
memleak.exe : memleak.cpp memleak_func.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@



####### End of Makefile #######
