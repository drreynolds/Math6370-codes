###################################################################
#  Makefile for Example: dot-product
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 4370/6370
#  18 October 2021
###################################################################

# compilers/flags
CXX = g++
CXXFLAGS = -O2 -std=c++11 -fopenmp
ifndef RAJA
  RAJA = /usr/local/raja-0.12.1
endif

# makefile targets
all : dot_prod_cpp.exe

dot_prod_cpp.exe : dot_prod.cpp get_time.c
	$(CXX) $(CXXFLAGS) -I$(RAJA)/include -L$(RAJA)/lib -lRAJA dot_prod.cpp get_time.c -o $@

clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
