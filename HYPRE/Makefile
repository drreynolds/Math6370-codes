###################################################################
#  Makefile for hypre test code
#  Daniel R. Reynolds
###################################################################

# C compiler
CXX = mpicxx -O0 -g

# HYPRE installation location
ifeq ($(HOST), asymptote)
    HYPREROOT  = /opt/local/hypre-2.24.0
else
    HYPREROOT  = /usr/local/hypre-2.20.0
endif
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
	\rm -rf *.exe *.dSYM

##### End of Makefile #####
