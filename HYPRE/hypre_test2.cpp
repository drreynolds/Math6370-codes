/*----------------------------------------------------------------
 Daniel R. Reynolds @ SMU Mathematics
 -----------------------------------------------------------------
 Copyright 2017, All rights reserved
 -----------------------------------------------------------------
 Description: 
   This example sets up (and solves) a scalar-valued Poisson-like 
   problem, with homogeneous Dirichlet boundary conditions.  
 
   Specifically, we set up and solve the linear system
           (I + L)v = w,
   where L is a standard 2D 5-pt Laplace stencil operator
   (anisotropic diffusion coefficients), and w is a smooth RHS 
   vector.  The Laplace operator is discretized using a simple 
   second-order centered difference approximation on a 
   cell-centered finite-volume grid.

   This routine uses HYPRE's BiCGSTAB solver, along with PFMG as 
   the preconditioner.  We use BiCGSTAB because the Dirichlet 
   boundary conditions cause trouble for the PCG solver (we 
   still use the input parameter PCGmaxit).
 -------------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "HYPRE_struct_ls.h"


// main program
int main(int argc, char *argv[]) {

  // local/reused variables
  int ierr, i, j;

  // Initialize MPI environment
  ierr = MPI_Init(&argc, &argv);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"Error in MPI_Init!\n");
    return 1;
  }

  int numprocs;
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"Error in MPI_Comm_size!\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int myid;
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"Error in MPI_Comm_rank!\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // set output processor flag
  bool outproc = (myid == 0);

  // Start main program timer
  double stime1 = MPI_Wtime();

  // Start initialization timer
  double stime2 = MPI_Wtime();

  // Input test parameters from file
  int Nx, Ny, Px, Py, BiCGSTABmaxit, PFMGmaxit, relch, rlxtype, npre, npost;
  double mu_x, mu_y, delta, xL, xR, yL, yR;
  if (myid == 0) {
    FILE* FID = fopen("input_params2.txt","r");
    ierr = fscanf(FID,"Nx = %i\n", &Nx);
    ierr = fscanf(FID,"Ny = %i\n", &Ny);
    ierr = fscanf(FID,"Px = %i\n", &Px);
    ierr = fscanf(FID,"Py = %i\n", &Py);
    ierr = fscanf(FID,"mu_x = %lf\n", &mu_x);
    ierr = fscanf(FID,"mu_y = %lf\n", &mu_y);
    ierr = fscanf(FID,"xL = %lf\n", &xL);
    ierr = fscanf(FID,"xR = %lf\n", &xR);
    ierr = fscanf(FID,"yL = %lf\n", &yL);
    ierr = fscanf(FID,"yR = %lf\n", &yR);
    ierr = fscanf(FID,"tol = %lf\n", &delta);
    ierr = fscanf(FID,"PCGmaxit = %i\n", &BiCGSTABmaxit);
    ierr = fscanf(FID,"PFMGmaxit = %i\n", &PFMGmaxit);
    ierr = fscanf(FID,"relch = %i\n", &relch);
    ierr = fscanf(FID,"rlxtype = %i\n", &rlxtype);
    ierr = fscanf(FID,"npre = %i\n", &npre);
    ierr = fscanf(FID,"npost = %i\n", &npost);
    fclose(FID);
  }
  
  // Broadcast inputs to all procs
  int ibuff[] = {Nx, Ny, Px, Py, BiCGSTABmaxit, PFMGmaxit, relch, rlxtype, npre, npost};
  double dbuff[] = {mu_x, mu_y, delta, xL, xR, yL, yR};
  ierr = MPI_Bcast(ibuff, 10, MPI_INT, 0, MPI_COMM_WORLD);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"Error in MPI_Bcast!\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  ierr = MPI_Bcast(dbuff, 7, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"Error in MPI_Bcast!\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  Nx = ibuff[0];
  Ny = ibuff[1];
  Px = ibuff[2];
  Py = ibuff[3];
  BiCGSTABmaxit = ibuff[4];
  PFMGmaxit = ibuff[5];
  relch = ibuff[6];
  rlxtype = ibuff[7];
  npre = ibuff[8];
  npost = ibuff[9];
  mu_x = dbuff[0];
  mu_y = dbuff[1];
  delta = dbuff[2];
  xL = dbuff[3];
  xR = dbuff[4];
  yL = dbuff[5];
  yR = dbuff[6];

  // Output test information to stdout
  if (outproc) {
    printf("\n hypre_test problem parameters:\n");
    printf("   Nx = %i\n", Nx);
    printf("   Ny = %i\n", Ny);
    printf("   Px = %i\n", Px);
    printf("   Py = %i\n", Py);
    printf("   mu_x = %g\n", mu_x);
    printf("   mu_y = %g\n", mu_y);
    printf("   xL = %g\n", xL);
    printf("   xR = %g\n", xR);
    printf("   yL = %g\n", yL);
    printf("   yR = %g\n", yR);
    printf("   tol = %g\n", delta);
    printf("   PCGmaxit = %i\n", BiCGSTABmaxit);
    printf("   PFMGmaxit = %i\n", PFMGmaxit);
    printf("   relch = %i\n", relch);
    printf("   rlxtype = %i\n", rlxtype);
    printf("   npre = %i\n", npre);
    printf("   npost = %i\n\n", npost);
  }

  // check for legal Px, Py inputs
  if (Px*Py != numprocs) {
    fprintf(stderr,"hypre_test error: illegal parallel decomposition!\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // check for a legal 2D domain
  if ((xL >= xR) || (yL >= yR)) {
    fprintf(stderr,"hypre_test error: illegal domain!\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // set up 2D Cartesian communicator
  int pdims[] = {Px, Py};
  int periods[] = {0, 0};
  MPI_Comm comm;
  ierr = MPI_Cart_create(MPI_COMM_WORLD, 2, pdims, periods, 0, &comm);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Cart_create = %i\n",ierr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  int pcoords[2];
  ierr = MPI_Cart_coords(comm, myid, 2, pcoords);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Cart_coords = %i\n",ierr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // determine local extents
  int ilower[] = {pcoords[0]*Nx/Px, pcoords[1]*Ny/Py};
  int iupper[] = {(pcoords[0]+1)*Nx/Px-1, (pcoords[1]+1)*Ny/Py-1};

  // set grid information
  int NxLoc = iupper[0]-ilower[0]+1;
  int NyLoc = iupper[1]-ilower[1]+1;
  double dx = (xR-xL)/Nx;
  double dy = (yR-yL)/Ny;
  double xLLoc = xL + ilower[0]*dx;
  double yLLoc = yL + ilower[1]*dy;
  double xRLoc = xLLoc + NxLoc*dx;
  double yRLoc = yLLoc + NyLoc*dy;

  // all procs output relevant subdomain information (for debugging)
  printf(" p%i [coords = (%i,%i)]:  subdomain [%g,%g]x[%g,%g],  extents [%i,%i]x[%i,%i]\n",
	 myid, pcoords[0], pcoords[1], xLLoc, xRLoc, yLLoc, yRLoc, 
	 ilower[0], iupper[0], ilower[1], iupper[1]);

  // Allocate vector memory, and set initial guess values 
  double *xx = new double[NxLoc*NyLoc];
  double *bb = new double[NxLoc*NyLoc];
  if ((xx==NULL) || (bb==NULL)) {
    fprintf(stderr,"hypre_test error: failed to allocate array data\n");
    MPI_Abort(comm, 1);
  }
  for (i=0; i<NxLoc*NyLoc; i++)  xx[i] = 0.0;
  
  // set rhs vector
  double xloc, yloc;
  double pi=3.14159265358979323846;
  for (j=0; j<NyLoc; j++) {
    yloc = yLLoc + (j+0.5)*dy;
    for (i=0; i<NxLoc; i++) {
      xloc = xLLoc + (i+0.5)*dx;
      bb[j*NxLoc + i] = sin(2.0*xloc*pi/xR) + sin(2.0*yloc*pi/yR);
    }
  }

  // HYPRE initialization

  // set up the grid
  //   create grid object
  HYPRE_StructGrid grid;
  HYPRE_StructGridCreate(comm, 2, &grid);

  //   set my grid extents as if we have one part with multiple boxes:
  //     each processor describes it's own global extents
  HYPRE_StructGridSetExtents(grid, ilower, iupper);

  //   assemble the grid
  HYPRE_StructGridAssemble(grid);
  
  //   set up the stencil
  HYPRE_StructStencil stencil;
  HYPRE_StructStencilCreate(2, 5, &stencil);

  //   set stencil entries
  int offset[2];
  //     dependency to bottom
  offset[0] = 0;  offset[1] = -1;
  HYPRE_StructStencilSetElement(stencil, 0, offset);
  //     dependency to left
  offset[0] = -1;  offset[1] = 0;
  HYPRE_StructStencilSetElement(stencil, 1, offset);
  //     dependency to self
  offset[0] = 0;  offset[1] = 0;
  HYPRE_StructStencilSetElement(stencil, 2, offset);
  //     dependency to right
  offset[0] = 1;  offset[1] = 0;
  HYPRE_StructStencilSetElement(stencil, 3, offset);
  //     dependency to top
  offset[0] = 0;  offset[1] = 1;
  HYPRE_StructStencilSetElement(stencil, 4, offset);
 

  // set up the matrix
  //   create Matrix object
  HYPRE_StructMatrix  A;
  HYPRE_StructMatrixCreate(comm, grid, stencil, &A);

  //   initialize matrix
  HYPRE_StructMatrixInitialize(A);

  //   set matrix values over grid
  //     vals holds template for (I-Laplace) stencil
  double vals[5] = {-mu_y/dy/dy, -mu_x/dx/dx, 
  		     1.0 + 2.0*(mu_x/dx/dx + mu_y/dy/dy),
		    -mu_x/dx/dx, -mu_y/dy/dy};
  int entries[5] = {0, 1, 2, 3, 4};

  //    loop over domain interior, filling in entries
  int index[2];
  int ix, iy, is, ie, js, je;
  is = (ilower[0] == 0) ? 1 : 0;
  ie = (iupper[0] == Nx-1) ? NxLoc-1 : NxLoc;
  js = (ilower[1] == 0) ? 1 : 0;
  je = (iupper[1] == Ny-1) ? NyLoc-1 : NyLoc;
  for (iy=js; iy<je; iy++) {
    index[1] = ilower[1] + iy;
    for (ix=is; ix<ie; ix++) {
      index[0] = ilower[0] + ix;
      HYPRE_StructMatrixSetValues(A, index, 5, entries, vals);
    } 
  }

  //    loop over boundaries, setting Dirichlet condition
  if (ilower[1] == 0) {
    iy=0;         // bottom
    index[1] = ilower[1] + iy;
    double Bvals[4] = {-mu_x/dx/dx, 
                       1.0 + 2.0*mu_x/dx/dx + 3.0*mu_y/dy/dy,
                       -mu_x/dx/dx, -mu_y/dy/dy};
    int Bentries[4] = {1, 2, 3, 4};
    for (ix=0; ix<NxLoc; ix++) {
      index[0] = ilower[0] + ix;
      HYPRE_StructMatrixSetValues(A, index, 4, Bentries, Bvals);
    }
  }
  if (ilower[0] == 0) {
    ix=0;         // left
    index[0] = ilower[0] + ix;
    double Bvals[4] = {-mu_y/dy/dy, 
                       1.0 + 3.0*mu_x/dx/dx + 2.0*mu_y/dy/dy,
                       -mu_x/dx/dx, -mu_y/dy/dy};
    int Bentries[4] = {0, 2, 3, 4};
    for (iy=0; iy<NyLoc; iy++) {
      index[1] = ilower[1] + iy;
      HYPRE_StructMatrixSetValues(A, index, 4, Bentries, Bvals);
    }
  }
  if (iupper[0] == Nx-1) {
    ix=NxLoc-1;   // right
    index[0] = ilower[0] + ix;
    double Bvals[4] = {-mu_y/dy/dy, -mu_x/dx/dx, 
                       1.0 + 3.0*mu_x/dx/dx + 2.0*mu_y/dy/dy,
                       -mu_y/dy/dy};
    int Bentries[4] = {0, 1, 2, 4};
    for (iy=0; iy<NyLoc; iy++) {
      index[1] = ilower[1] + iy;
      HYPRE_StructMatrixSetValues(A, index, 4, Bentries, Bvals);
    }
  }
  if (iupper[1] == Ny-1) {
    iy=NyLoc-1;   // top
    index[1] = ilower[1] + iy;
    double Bvals[4] = {-mu_y/dy/dy, -mu_x/dx/dx, 
                       1.0 + 2.0*mu_x/dx/dx + 2.0*mu_y/dy/dy,
                       -mu_x/dx/dx};
    int Bentries[4] = {0, 1, 2, 3};
    for (ix=0; ix<NxLoc; ix++) {
      index[0] = ilower[0] + ix;
      HYPRE_StructMatrixSetValues(A, index, 4, Bentries, Bvals);
    }
  }

  //     assemble matrix
  HYPRE_StructMatrixAssemble(A);
    
  // Stop initialization timer
  double ftime2 = MPI_Wtime();

  // Start overall solver timer
  double stime3 = MPI_Wtime();

  // create the Struct vectors
  HYPRE_StructVector bvec, xvec;
  HYPRE_StructVectorCreate(comm, grid, &bvec);
  HYPRE_StructVectorCreate(comm, grid, &xvec);

  // initialize vectors
  HYPRE_StructVectorInitialize(bvec);
  HYPRE_StructVectorInitialize(xvec);  

  // convert rhs, solution vectors to HYPRE format     
  //    insert rhs vector entries into HYPRE vector b  
  //    and sol vector guesses into HYPRE vector x     
  for (iy=0; iy<NyLoc; iy++) {
    index[1] = ilower[1] + iy;
    for (ix=0; ix<NxLoc; ix++) {
      index[0] = ilower[0] + ix;
      HYPRE_StructVectorSetValues(bvec, index, bb[iy*NxLoc + ix]);
      HYPRE_StructVectorSetValues(xvec, index, xx[iy*NxLoc + ix]);
    }
  }

  //    assemble vectors
  HYPRE_StructVectorAssemble(xvec);
  HYPRE_StructVectorAssemble(bvec);

  // set up the preconditioner
  HYPRE_StructSolver precond;
  HYPRE_StructPFMGCreate(comm, &precond);
  HYPRE_StructPFMGSetMaxIter(precond, PFMGmaxit);
  HYPRE_StructPFMGSetRelaxType(precond, rlxtype);
  HYPRE_StructPFMGSetNumPreRelax(precond, npre);
  HYPRE_StructPFMGSetNumPostRelax(precond, npost);

  // set up the solver
  HYPRE_StructSolver solver;
  HYPRE_StructBiCGSTABCreate(comm, &solver);
  HYPRE_StructBiCGSTABSetMaxIter(solver, BiCGSTABmaxit);
  if (relch)
    HYPRE_StructBiCGSTABSetTol(solver, delta);
  else
    HYPRE_StructBiCGSTABSetAbsoluteTol(solver, delta);
  HYPRE_StructBiCGSTABSetLogging(solver, 1);
  HYPRE_StructBiCGSTABSetPrecond(solver,
			    (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSolve,
			    (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSetup,
			    precond);
  HYPRE_StructBiCGSTABSetup(solver, A, bvec, xvec);

  // Solve system, get iteration counts
  if (outproc)
    printf("\n hypre_test: solving matrix system \n");

  // Start hypre timer
  double stime4 = MPI_Wtime();

  // solve the linear system
  HYPRE_StructBiCGSTABSolve(solver, A, bvec, xvec);

  // Stop hypre timer
  double ftime4 = MPI_Wtime();

  // extract solver statistics
  double finalresid;
  HYPRE_StructBiCGSTABGetFinalRelativeResidualNorm(solver, &finalresid);
  int BiCGSTABits, PFMGits;
  HYPRE_StructBiCGSTABGetNumIterations(solver, &BiCGSTABits);
  HYPRE_StructPFMGGetNumIterations(precond, &PFMGits);
  if (outproc)
    printf("   lin resid = %.1e (tol = %.1e), BiCGSTAB its = %i, PFMG its = %i\n",
	   finalresid, delta, BiCGSTABits, PFMGits);

  // extract solution values
  double val;
  for (iy=0; iy<NyLoc; iy++) {
    index[1] = ilower[1] + iy;
    for (ix=0; ix<NxLoc; ix++) {
      index[0] = ilower[0] + ix;
      HYPRE_StructVectorGetValues(xvec, index, &val);
      xx[iy*NxLoc + ix] = val;
    }
  }

  // Stop overall timer
  double ftime1 = MPI_Wtime();

  // Output runtimes 
  if (outproc) {
    printf("\n hypre_test runtimes:\n");
    printf("   initialization time = %g\n",ftime2-stime2);
    printf("   overall solve time  = %g\n",ftime1-stime3);
    printf("   BiCGSTAB solve time = %g\n",ftime4-stime4);
    printf("   overall run time    = %g\n\n",ftime1-stime1);
  }

  // destroy HYPRE solver structures
  HYPRE_StructBiCGSTABDestroy(solver);
  HYPRE_StructPFMGDestroy(precond);
  HYPRE_StructVectorDestroy(bvec);
  HYPRE_StructVectorDestroy(xvec);
  HYPRE_StructMatrixDestroy(A);
  HYPRE_StructGridDestroy(grid);
  HYPRE_StructStencilDestroy(stencil);

  // free temporary vector memory
  delete[] xx;
  delete[] bb;

  // finalize MPI
  MPI_Finalize();

  return(0);
}

