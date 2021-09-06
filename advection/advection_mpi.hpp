/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   11 May 2017 */

#ifndef _ADVECTION_MPI_HPP
#define _ADVECTION_MPI_HPP

#include "mpi.h"

// simple macro to map a 2D index to a 1D address space
#define idx(i,j,nx)   ((j)*(nx)+(i))


// Parallelism class
class parallel_decomp {

 public:

  MPI_Comm comm;
  int pdims[2];
  int periodic[2];
  int pcoords[2];
  int nxloc;
  int nyloc;
  int myid;
  int numprocs;
  int nbE;
  int nbW;
  int nbN;
  int nbS;
  double *v2recvE;
  double *v3recvN;
  double *v2sendW;
  double *v3sendS;
  double *v1recvW;
  double *v1recvS;
  double *v1sendE;
  double *v1sendN;

  parallel_decomp() {
    nbE = nbW = nbN = nbS = -1;
    v2recvE = v3recvN = v2sendW = v3sendS = v1recvW = v1recvS = v1sendE = v1sendN = NULL;
  }
  ~parallel_decomp() {
    if (v2recvE != NULL)  delete[] v2recvE;
    if (v3recvN != NULL)  delete[] v3recvN;
    if (v2sendW != NULL)  delete[] v2sendW;
    if (v3sendS != NULL)  delete[] v3sendS;
    if (v1recvW != NULL)  delete[] v1recvW;
    if (v1recvS != NULL)  delete[] v1recvS;
    if (v1sendE != NULL)  delete[] v1sendE;
    if (v1sendN != NULL)  delete[] v1sendN;
  }
  int Communication1(double *v2, double *v3);
  int Communication2(double *v1);
  
};  // end parallel_decomp



// Prototypes for utility routines
void initialize(double* u, double* v1, double* v2, double* v3, double c, 
		double dx, double dy, int is, int ie, int js, int je);

void output(double* u, double t, int nx, int ny, int noutput, parallel_decomp& p2d);

void check_err(int ierr, MPI_Comm comm, const char* fname);

#endif
