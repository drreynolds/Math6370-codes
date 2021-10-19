###################################################################
#  Makefile for Example: matrix-vector product
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
all : matvec_row_cpp.exe matvec_col_cpp.exe

matvec_row_cpp.exe : matvec_row.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

matvec_col_cpp.exe : matvec_col.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@


clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~ b*.txt


####### End of Makefile #######
