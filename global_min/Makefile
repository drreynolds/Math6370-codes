###################################################################
#  Makefile for Example: global minimization
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 4370 / 6370
###################################################################

# set compilers based on machine
CXX = g++
CC  = cc
F90 = gfortran

# flags
CXXFLAGS = -O2
CFLAGS   = -O2
FFLAGS   = -O2


# makefile targets
all : all_cpp all_c all_f90

all_cpp : glob_min_cpp.exe glob_min2_cpp.exe

all_c : glob_min_c.exe glob_min2_c.exe

all_f90 : glob_min_f90.exe glob_min2_f90.exe


glob_min_cpp.exe : glob_min.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

glob_min2_cpp.exe : glob_min2.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@


glob_min_c.exe : glob_min.c get_time.c
	$(CC) $(CFLAGS) $^ -lm -o $@

glob_min2_c.exe : glob_min2.c get_time.c
	$(CC) $(CFLAGS) $^ -lm -o $@


glob_min_f90.exe : glob_min.f90 get_time.f90
	$(F90) $(FFLAGS) $^ -o $@

glob_min2_f90.exe : glob_min2.f90 get_time.f90
	$(F90) $(FFLAGS) $^ -o $@


clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
