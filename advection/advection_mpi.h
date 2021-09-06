/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   11 May 2017 */

#ifndef _ADVECTION_H
#define _ADVECTION_H

#include "mpi.h"

// simple macro to map a 2D index to a 1D address space
#define idx(i,j,nx)   ((j)*(nx)+(i))


// Parallelism structure
typedef struct _parallel_decomp {

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

} parallel_decomp;


// Prototypes
void create_parallel_decomp(parallel_decomp *p2d);

void free_parallel_decomp(parallel_decomp *p2d);

int Communication1(double *v2, double *v3, parallel_decomp *p2d);

int Communication2(double *v1, parallel_decomp *p2d);

void initialize(double *u, double *v1, double *v2, double *v3, double c, 
                double dx, double dy, int is, int ie, int js, int je);

void output(double *u, double t, int nx, int ny, int noutput, parallel_decomp *p2d);

void check_err(int ierr, MPI_Comm comm, const char* fname);

#endif
