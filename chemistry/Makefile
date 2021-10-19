###################################################################
#  Makefile for Example: computing chemical equlibria
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
all : chemistry_cpp.exe

chemistry_cpp.exe : chemistry.cpp chem_solver.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
