###################################################################
#  Makefile for Example: computing chemical equlibria
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 4370/6370
#  7 February 2015
###################################################################

# F90
F90 = gfortran
MPIF90 = mpif90
FFLAGS = -O2

# C++
CXX = g++
MPICXX = mpicxx
CXXFLAGS = -O2

# C
CC = gcc
MPICC = mpicc
CFLAGS = -O2

# OpenMP
OMPFLAG = -fopenmp


# makefile targets
all : cpp c f90


cpp : chemistry_cpp.exe chemistry_mpi_cpp.exe \
      chemistry_mpi2_cpp.exe chemistry_hybrid_cpp.exe \
      #chemistry_mpi3_cpp.exe 

chemistry_cpp.exe : chemistry.cpp chem_solver.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

chemistry_mpi_cpp.exe : chemistry_mpi.cpp chem_solver.cpp
	$(MPICXX) $(CXXFLAGS) $^ -lm -o $@

chemistry_mpi2_cpp.exe : chemistry_mpi2.cpp chem_solver.cpp
	$(MPICXX) $(CXXFLAGS) $^ -lm -o $@

chemistry_mpi3_cpp.exe : chemistry_mpi3.cpp chem_solver.cpp 
	$(MPICXX) $(CXXFLAGS) $^ -lm -o $@

chemistry_hybrid_cpp.exe : chemistry_hybrid.cpp chem_solver.cpp 
	$(MPICXX) $(CXXFLAGS) $(OMPFLAG) $^ -o $@



c : chemistry_c.exe chemistry_mpi_c.exe \
    chemistry_mpi2_c.exe #chemistry_mpi3_c.exe 

chemistry_c.exe : chemistry.c chem_solver.c get_time.c
	$(CC) $(CFLAGS) $^ -lm -o $@

chemistry_mpi_c.exe : chemistry_mpi.c chem_solver.c
	$(MPICC) $(CFLAGS) $^ -lm -o $@

chemistry_mpi2_c.exe : chemistry_mpi2.c chem_solver.c
	$(MPICC) $(CFLAGS) $^ -lm -o $@

chemistry_mpi3_c.exe : chemistry_mpi3.c chem_solver.c
	$(MPICC) $(CFLAGS) $^ -lm -o $@


f90 : chemistry_f90.exe chemistry_mpi_f90.exe \
      chemistry_mpi2_f90.exe chemistry_mpi3_f90.exe

chemistry_f90.exe : chemistry.f90 chem_solver.f90 get_time.f90
	$(F90) $(FFLAGS) $^ -o $@

chemistry_mpi_f90.exe : chemistry_mpi.f90 chem_solver.f90
	$(MPIF90) $(FFLAGS) $^ -o $@

chemistry_mpi2_f90.exe : chemistry_mpi2.f90 chem_solver.f90
	$(MPIF90) $(FFLAGS) $^ -o $@

chemistry_mpi3_f90.exe : chemistry_mpi3.f90 chem_solver.f90
	$(MPIF90) $(FFLAGS) $^ -o $@


clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
