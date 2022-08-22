###################################################################
#  Makefile for Example: computing chemical equlibria
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 4370/6370
#  7 February 2015
###################################################################

# set compilers based on machine
ifeq ($(HOST), asymptote)
  CXX = g++-12
  CC  = gcc-12
  F90 = gfortran-12
else
  CXX = g++
  CC  = gcc
  F90 = gfortran
endif

# flags
CXXFLAGS = -O2
CFLAGS   = -O2
FFLAGS   = -O2

OMPFLAGS = -fopenmp
AUTOPAR  = -floop-parallelize-all -ftree-parallelize-loops=4


# makefile targets
all : all_cpp all_c all_f90

all_cpp : chemistry_cpp.exe chemistry_auto_cpp.exe \
      chemistry_omp_cpp.exe chemistry_omp_orphan_cpp.exe

all_c : chemistry_c.exe chemistry_auto_c.exe \
      chemistry_omp_c.exe chemistry_omp_orphan_c.exe

all_f90 : chemistry_f90.exe chemistry_auto_f90.exe \
      chemistry_omp_f90.exe chemistry_omp_orphan_f90.exe


chemistry_cpp.exe : chemistry.cpp chem_solver.cpp get_time.c
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

chemistry_auto_cpp.exe : chemistry.cpp chem_solver.cpp get_time.c
	$(CXX) $(CXXFLAGS) $(AUTOPAR) $^ -lm -o $@

chemistry_omp_cpp.exe : chemistry_omp.cpp chem_solver.cpp get_time.c
	$(CXX) $(OMPFLAGS) $(CXXFLAGS) $^ -lm -o $@

chemistry_omp_orphan_cpp.exe : chemistry_omp2.cpp chem_solver_orphan.cpp get_time.c
	$(CXX) $(OMPFLAGS) $(CXXFLAGS) $^ -lm -o $@


chemistry_c.exe : chemistry.c chem_solver.c get_time.c
	$(CC) $(CFLAGS) $^ -lm -o $@

chemistry_auto_c.exe : chemistry.c chem_solver.c get_time.c
	$(CC) $(CFLAGS) $(AUTOPAR) $^ -lm -o $@

chemistry_omp_c.exe : chemistry_omp.c chem_solver.c get_time.c
	$(CC) $(OMPFLAGS) $(CFLAGS) $^ -lm -o $@

chemistry_omp_orphan_c.exe : chemistry_omp2.c chem_solver_orphan.c get_time.c
	$(CC) $(OMPFLAGS) $(CCFLAGS) $^ -lm -o $@


chemistry_f90.exe : chemistry.f90 chem_solver.f90 get_time.f90
	$(F90) $(FFLAGS) $^ -o $@

chemistry_auto_f90.exe : chemistry.f90 chem_solver.f90 get_time.f90
	$(F90) $(FFLAGS) $(AUTOPAR) $^ -o $@

chemistry_omp_f90.exe : chemistry_omp.f90 chem_solver.f90 get_time.f90
	$(F90) $(OMPFLAGS) $(FFLAGS) $^ -o $@

chemistry_omp_orphan_f90.exe : chemistry_omp2.f90 chem_solver_orphan.f90 get_time.f90
	$(F90) $(OMPFLAGS) $(FFLAGS) $^ -o $@


clean :
	\rm -f *.o *.out

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
