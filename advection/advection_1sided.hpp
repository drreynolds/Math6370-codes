/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

#ifndef _ADVECTION_1SIDED_HPP
#define _ADVECTION_1SIDED_HPP

#include <iostream>
#include "mpi.h"

// simple macro to map a 2D index to a 1D address space
#define idx(i,j,nx)   ((j)*(nx)+(i))

// Parallelism class
class parallel_decomp {

 public:

  MPI_Comm comm;
  MPI_Win win_v1W;
  MPI_Win win_v1S;
  MPI_Win win_v2E;
  MPI_Win win_v3N;
  MPI_Group groupE;
  MPI_Group groupW;
  MPI_Group groupN;
  MPI_Group groupS;
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
  double *v2_recvE;
  double *v3_recvN;
  double *v2_sendW;
  double *v3_sendS;
  double *v1_recvW;
  double *v1_recvS;
  double *v1_sendE;
  double *v1_sendN;

  parallel_decomp() {
    comm = MPI_COMM_NULL;
    win_v1W = win_v1S = win_v2E = win_v3N = MPI_WIN_NULL;
    groupE = groupW = groupN = groupS = MPI_GROUP_NULL;
    pdims[0] = 0;
    pdims[1] = 0;
    periodic[0] = 0;
    periodic[1] = 0;
    pcoords[0] = 0;
    pcoords[1] = 0;
    is = ie = js = je = nxloc = nyloc = myid = numprocs = -1;
    nbE = nbW = nbN = nbS = MPI_PROC_NULL;
    v2_recvE = v3_recvN = v2_sendW = v3_sendS = v1_recvW =
      v1_recvS = v1_sendE = v1_sendN = nullptr;
  }
  void free() {
    if (this->v2_recvE != nullptr) {
      delete[] this->v2_recvE;
      this->v2_recvE = nullptr;
    }
    if (this->v3_recvN != nullptr) {
      delete[] this->v3_recvN;
      this->v3_recvN = nullptr;
    }
    if (this->v2_sendW != nullptr) {
      delete[] this->v2_sendW;
      this->v2_sendW = nullptr;
    }
    if (this->v3_sendS != nullptr) {
      delete[] this->v3_sendS;
      this->v3_sendS = nullptr;
    }
    if (this->v1_recvW != nullptr) {
      delete[] this->v1_recvW;
      this->v1_recvW = nullptr;
    }
    if (this->v1_recvS != nullptr) {
      delete[] this->v1_recvS;
      this->v1_recvS = nullptr;
    }
    if (this->v1_sendE != nullptr) {
      delete[] this->v1_sendE;
      this->v1_sendE = nullptr;
    }
    if (this->v1_sendN != nullptr) {
      delete[] this->v1_sendN;
      this->v1_sendN = nullptr;
    }
    if (this->win_v1W != MPI_WIN_NULL) {
      MPI_Win_free(&(this->win_v1W));
      this->win_v1W = MPI_WIN_NULL;
    }
    if (this->win_v1S != MPI_WIN_NULL) {
      MPI_Win_free(&(this->win_v1S));
      this->win_v1S = MPI_WIN_NULL;
    }
    if (this->win_v2E != MPI_WIN_NULL) {
      MPI_Win_free(&(this->win_v2E));
      this->win_v2E = MPI_WIN_NULL;
    }
    if (this->win_v3N != MPI_WIN_NULL) {
      MPI_Win_free(&(this->win_v3N));
      this->win_v3N = MPI_WIN_NULL;
    }
    if (this->groupE != MPI_GROUP_NULL) {
      MPI_Group_free(&(this->groupE));
      this->groupE = MPI_GROUP_NULL;
    }
    if (this->groupW != MPI_GROUP_NULL) {
      MPI_Group_free(&(this->groupW));
      this->groupW = MPI_GROUP_NULL;
    }
    if (this->groupN != MPI_GROUP_NULL) {
      MPI_Group_free(&(this->groupN));
      this->groupN = MPI_GROUP_NULL;
    }
    if (this->groupS != MPI_GROUP_NULL) {
      MPI_Group_free(&(this->groupS));
      this->groupS = MPI_GROUP_NULL;
    }
  }
  ~parallel_decomp() {
    this->free();
  }
  void setup(int nx, int ny);
  void start_v2v3();
  void put_v2v3(double* v2, double* v3);
  void finish_v2v3();
  void start_v1();
  void put_v1(double* v1);
  void finish_v1();

};  // end parallel_decomp


// Prototypes for utility routines
void initialize(double* u, double* v1, double* v2, double* v3,
		double c, double dx, double dy, parallel_decomp& p2d);

void output(double* u, double t, int nx, int ny, int noutput, parallel_decomp& p2d);

// error checking routine for successful MPI calls
void check_err(const int ierr, MPI_Comm comm, const char* fname);

#endif
