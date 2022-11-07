/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

#ifndef _ADVECTION_MPI_HPP
#define _ADVECTION_MPI_HPP

#include <iostream>
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
  int is;
  int ie;
  int js;
  int je;
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
    comm = MPI_COMM_NULL;
    pdims[0] = 0;
    pdims[1] = 0;
    periodic[0] = 0;
    periodic[1] = 0;
    pcoords[0] = 0;
    pcoords[1] = 0;
    is = ie = js = je = nxloc = nyloc = myid = numprocs = -1;
    nbE = nbW = nbN = nbS = MPI_PROC_NULL;
    v2recvE = v3recvN = v2sendW = v3sendS = v1recvW = v1recvS = v1sendE = v1sendN = nullptr;
  }
  void free() {
    if (this->v2recvE != nullptr) {
      delete[] this->v2recvE;
      this->v2recvE = nullptr;
    }
    if (this->v3recvN != nullptr) {
      delete[] this->v3recvN;
      this->v3recvN = nullptr;
    }
    if (this->v2sendW != nullptr) {
      delete[] this->v2sendW;
      this->v2sendW = nullptr;
    }
    if (this->v3sendS != nullptr) {
      delete[] this->v3sendS;
      this->v3sendS = nullptr;
    }
    if (this->v1recvW != nullptr) {
      delete[] this->v1recvW;
      this->v1recvW = nullptr;
    }
    if (this->v1recvS != nullptr) {
      delete[] this->v1recvS;
      this->v1recvS = nullptr;
    }
    if (this->v1sendE != nullptr) {
      delete[] this->v1sendE;
      this->v1sendE = nullptr;
    }
    if (this->v1sendN != nullptr) {
      delete[] this->v1sendN;
      this->v1sendN = nullptr;
    }
  }
  ~parallel_decomp() {
    this->free();
  }
  void setup(int nx, int ny);
  int Communication1(double *v2, double *v3);
  int Communication2(double *v1);

};  // end parallel_decomp


// Prototypes for utility routines
void initialize(double* u, double* v1, double* v2, double* v3,
                double c, double dx, double dy, parallel_decomp& p2d);

void output(double* u, double t, int nx, int ny, int noutput, parallel_decomp& p2d);

// error checking routine for successful MPI calls
void check_err(const int ierr, MPI_Comm comm, const char* fname);

#endif
