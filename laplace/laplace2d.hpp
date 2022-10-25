/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

#ifndef _LAPLACE2D_HPP
#define _LAPLACE2D_HPP

// simple macro to map a 2D index to a 1D address space
#define idx(i,j,nx)   ((j)*(nx)+(i))

// Prototypes
int linresid2D(double *u, double *f, double *res, double &norm2,
	       int locN, int locM, double dx, double dy, MPI_Comm comm);

#endif
