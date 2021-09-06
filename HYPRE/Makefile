###################################################################
#  Makefile for hypre test code
#  Daniel R. Reynolds
###################################################################

# C compiler
CXX = mpicxx -O0 -g

# HYPRE installation location
HYPREROOT  = /usr/local/hypre-2.20.0
HYPRE_INCD = -I${HYPREROOT}/include
HYPRE_LIBD = -L${HYPREROOT}/lib
HYPRE_LIBS = -lHYPRE


##### Makefile targets #####

all : hypre_test.exe hypre_test2.exe hypre_test_fd.exe hypre_test2_fd.exe

hypre_test.exe : hypre_test.cpp
	${CXX} -o $@ $^ ${HYPRE_INCD} ${HYPRE_LIBD} ${HYPRE_LIBS}

hypre_test2.exe : hypre_test2.cpp
	${CXX} -o $@ $^ ${HYPRE_INCD} ${HYPRE_LIBD} ${HYPRE_LIBS}

hypre_test_fd.exe : hypre_test_fd.cpp
	${CXX} -o $@ $^ ${HYPRE_INCD} ${HYPRE_LIBD} ${HYPRE_LIBS}

hypre_test2_fd.exe : hypre_test2_fd.cpp
	${CXX} -o $@ $^ ${HYPRE_INCD} ${HYPRE_LIBD} ${HYPRE_LIBS}

clean:
	\rm -f *.o *~

realclean: clean
	\rm -f *.exe

##### End of Makefile #####
