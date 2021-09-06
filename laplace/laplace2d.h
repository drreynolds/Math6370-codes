/* Daniel R. Reynolds
   SMU Mathematics
   Math 6370
   28 February 2013 */

#ifndef _LAPLACE2D_H
#define _LAPLACE2D_H

// simple macro to map a 2D index to a 1D address space
#define idx(i,j,nx)   ((j)*(nx)+(i))

// Prototypes 
int linresid2D(double *u, double *f, double *res, double *norm2, 
	       int locN, int locM, double dx, double dy, MPI_Comm comm);

#endif
