###################################################################
#  Makefile for Example: matrix-vector product
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 6370
#  28 February 2013
###################################################################

MPIF90 = mpif90
FFLAGS = -O2

MPICXX = mpicxx
CFLAGS = -O2 -lm


# makefile targets
all : cpp f90


cpp : direct_cpp.exe iterative_cpp.exe

direct_cpp.exe: direct_main.cpp tridiag_solve.cpp linresid.cpp
	$(MPICXX) $(CFLAGS) -o $@ $^

iterative_cpp.exe: iterative_main.cpp jacobi_solve.cpp linresid.cpp
	$(MPICXX) $(CFLAGS) -o $@ $^ 


f90 : direct_f90.exe iterative_f90.exe

direct_f90.exe: direct_main.f90 tridiag_solve.f90 linresid.f90
	$(MPIF90) $(FFLAGS) -o $@ $^

iterative_f90.exe: iterative_main.f90 jacobi_solve.f90 linresid.f90
	$(MPIF90) $(FFLAGS) -o $@ $^


clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~ b*.txt


####### End of Makefile #######
