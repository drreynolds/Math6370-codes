###################################################################
#  Makefile for Example: HYPRE
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 4370 / 6370
###################################################################

# check for required environment variables
ifndef HYPRE_ROOT
    $(error Missing HYPRE_ROOT environment variable)
endif
ifndef MPICXX
    $(error Missing MPICXX environment variable)
endif

# set HYPRE-related variables
HYPRE_INCD = -I${HYPRE_ROOT}/include
HYPRE_LIBD = -L${HYPRE_ROOT}/lib -Wl,-rpath,${HYPRE_ROOT}/lib
HYPRE_LIBS = -lHYPRE


##### Makefile targets #####

all : hypre_test.exe hypre_test2.exe hypre_test_fd.exe hypre_test2_fd.exe

hypre_test.exe : hypre_test.cpp
	${MPICXX} -o $@ $^ ${HYPRE_INCD} ${HYPRE_LIBD} ${HYPRE_LIBS}

hypre_test2.exe : hypre_test2.cpp
	${MPICXX} -o $@ $^ ${HYPRE_INCD} ${HYPRE_LIBD} ${HYPRE_LIBS}

hypre_test_fd.exe : hypre_test_fd.cpp
	${MPICXX} -o $@ $^ ${HYPRE_INCD} ${HYPRE_LIBD} ${HYPRE_LIBS}

hypre_test2_fd.exe : hypre_test2_fd.cpp
	${MPICXX} -o $@ $^ ${HYPRE_INCD} ${HYPRE_LIBD} ${HYPRE_LIBS}

clean:
	\rm -f *.o *~

realclean: clean
	\rm -rf *.exe *.dSYM

##### End of Makefile #####
