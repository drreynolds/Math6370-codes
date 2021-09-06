###################################################################
#  Makefile for Example: Laplace residual
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 6370
#  28 February 2013
###################################################################

# F90
MPIF90 = mpif90
FFLAGS = -O2

# C++
MPICXX = mpicxx
CXXFLAGS = -O2 -lm

# C
MPICC = mpicc
CFLAGS = -O2


all : cpp c f90


cpp : driver_cpp.exe driver_async_cpp.exe

driver_cpp.exe: driver.cpp linresid2D_mpi_sync.cpp
	$(MPICXX) $(CXXFLAGS) -o $@ $^

driver_async_cpp.exe: driver.cpp linresid2D_mpi_async.cpp
	$(MPICXX) $(CXXFLAGS) -o $@ $^


c : driver_c.exe driver_async_c.exe

driver_c.exe: driver.c linresid2D_mpi_sync.c
	$(MPICC) $(CFLAGS) -o $@ $^ -lm

driver_async_c.exe: driver.c linresid2D_mpi_async.c
	$(MPICC) $(CFLAGS) -o $@ $^ -lm


f90 : driver_f90.exe driver_async_f90.exe

driver_f90.exe: driver.f90 linresid2D_mpi_sync.f90
	$(MPIF90) $(FFLAGS) -o $@ $^

driver_async_f90.exe: driver.f90 linresid2D_mpi_async.f90
	$(MPIF90) $(FFLAGS) -o $@ $^


clean :
	\rm -f *.o *.out resid.*

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
