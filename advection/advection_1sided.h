/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   30 April 2015 */

#include "mpi.h"
#include <stdio.h>

// simple macro to map a 2D index to a 1D address space
#define idx(i,j,nx)   ((j)*(nx)+(i))

// Prototypes
void initialize(double* u, double* v1, double* v2, double* v3, 
		double c, double dx, double dy, int is, 
		int js, int nxl, int nyl);

void output(double* u, double t, int nx, int ny, int nxl, 
	    int nyl, int noutput, MPI_Comm comm);

int setup_windows(int nxl, int nyl, double* v1_recvW, 
		  double* v1_recvS, double* v2_recvE, 
		  double* v3_recvN, MPI_Comm comm, 
		  MPI_Win& win_v1W, MPI_Win& win_v1S, 
		  MPI_Win& win_v2E, MPI_Win& win_v3N,
		  MPI_Group& groupW, MPI_Group& groupE, 
		  MPI_Group& groupS, MPI_Group& groupN);

int put_v2v3(double* v2, double* v3, double* v2_sendW, 
	     double* v3_sendS, int nxl, int nyl, MPI_Comm comm, 
	     MPI_Win& win_v2E, MPI_Win& win_v3N);

int put_v1(double* v1, double* v1_sendE, double* v1_sendN, int nxl, 
	   int nyl, MPI_Comm comm, MPI_Win& win_v1W, MPI_Win& win_v1S);

void check_err(int ierr, MPI_Comm comm, const char* fname);

