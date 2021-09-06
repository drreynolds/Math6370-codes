/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   11 May 2017 */

/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "advection_mpi.h"
#include "mpi.h"


/* Example routine to evolve the first-order 2D wave equations in time. */
int main(int argc, char* argv[]) {

  /* declarations */
  int nx, ny, nt, i, j, it, ierr, ibuff[3];
  double tstop, c, dtoutput;
  double v1_N, v1_S, v1_E, v1_W, v2_E, v2_W, v3_N, v3_S, rbuff[3];

  /* intialize MPI */
  parallel_decomp p2d;
  create_parallel_decomp(&p2d);
  
  ierr = MPI_Init(&argc, &argv);
  check_err(ierr, MPI_COMM_WORLD, "MPI_Init");

  ierr = MPI_Comm_size(MPI_COMM_WORLD, &(p2d.numprocs));
  check_err(ierr, MPI_COMM_WORLD, "MPI_Comm_size");

  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &(p2d.myid));
  check_err(ierr, MPI_COMM_WORLD, "MPI_Comm_size");

  double stime = MPI_Wtime();
  /* read problem parameters from input file (should be in this order):
        nx - number of nodes in x-direction
	ny - number of nodes in y-direction
	nt - number of time steps
	tstop - final time (will stop at nt or stop, whichever is 1st)
	c - wave speed
	dtoutput - time frequency for outputting solutions */
  if (p2d.myid == 0) {
    FILE *FID = fopen("input.txt","r");
    i = fscanf(FID," &inputs\n");
    i = fscanf(FID,"  nx = %i,\n", &nx);
    i = fscanf(FID,"  ny = %i,\n", &ny);
    i = fscanf(FID,"  nt = %i,\n", &nt);
    i = fscanf(FID,"  tstop = %lf,\n", &tstop);
    i = fscanf(FID,"  c = %lf,\n", &c);
    i = fscanf(FID,"  dtoutput = %lf,\n", &dtoutput);
    fclose(FID);

    printf("\nRunning wave problem:\n");
    printf("  nx = %i,  ny = %i\n", nx, ny);
    printf("  nt = %i,  tstop = %g\n", nt, tstop);
    printf("  c = %g\n", c);
    printf("  dtoutput = %g\n", dtoutput);
    printf("  nprocs = %i\n", p2d.numprocs);

    // fill buffers for other MPI processes
    ibuff[0] = nx;
    ibuff[1] = ny;
    ibuff[2] = nt;
    rbuff[0] = tstop;
    rbuff[1] = c;
    rbuff[2] = dtoutput;
  }
  /* broadcast parameters to all procs */
  ierr = MPI_Bcast(ibuff, 3, MPI_INT, 0, MPI_COMM_WORLD);
  check_err(ierr, MPI_COMM_WORLD, "MPI_Bcast");
  ierr = MPI_Bcast(rbuff, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  check_err(ierr, MPI_COMM_WORLD, "MPI_Bcast");
  nx = ibuff[0];
  ny = ibuff[1];
  nt = ibuff[2];
  tstop = rbuff[0];
  c = rbuff[1];
  dtoutput = rbuff[2];
  double ftime = MPI_Wtime();
  double iotime = ftime - stime;

  /* perform parallel decomposition */
  stime = MPI_Wtime();
  p2d.periodic[0] = 1;
  p2d.periodic[1] = 1;
  p2d.pdims[0] = 0;
  p2d.pdims[1] = 0;
  ierr = MPI_Dims_create(p2d.numprocs, 2, p2d.pdims);
  check_err(ierr, MPI_COMM_WORLD, "MPI_Dims_create");

  ierr = MPI_Cart_create(MPI_COMM_WORLD, 2, p2d.pdims, p2d.periodic, 0, &(p2d.comm));
  check_err(ierr, MPI_COMM_WORLD, "MPI_Cart_create");

  ierr = MPI_Cart_coords(p2d.comm, p2d.myid, 2, p2d.pcoords);
  check_err(ierr, p2d.comm, "MPI_Cart_get");

  int is = p2d.pcoords[0]*nx/p2d.pdims[0];
  int ie = (p2d.pcoords[0]+1)*nx/p2d.pdims[0]-1;
  int js = p2d.pcoords[1]*ny/p2d.pdims[1];
  int je = (p2d.pcoords[1]+1)*ny/p2d.pdims[1]-1;

  /* allocate arrays */
  p2d.nxloc = ie-is+1;
  p2d.nyloc = je-js+1;
  double *u  = (double *) malloc((p2d.nxloc*p2d.nyloc) * sizeof(double));
  double *v1 = (double *) malloc((p2d.nxloc*p2d.nyloc) * sizeof(double));
  double *v2 = (double *) malloc((p2d.nxloc*p2d.nyloc) * sizeof(double));
  double *v3 = (double *) malloc((p2d.nxloc*p2d.nyloc) * sizeof(double));

  /* set grid spacing */
  double dx = 1.0/nx;
  double dy = 1.0/ny;

  /* set initial conditions */
  initialize(u, v1, v2, v3, c, dx, dy, is, ie, js, je);
  ftime = MPI_Wtime();
  double inittime = ftime - stime;

  /* set initial time, output initial solution */
  double t = 0.0;
  double toutput = 0.0;
  int noutput = 0;
  if (p2d.myid == 0)
    printf("writing output file %i, step = %i, t = %g\n", noutput, 0, t);
  stime = MPI_Wtime();
  output(u, t, nx, ny, noutput, &p2d);
  ftime = MPI_Wtime();
  iotime += ftime - stime;

  /* set time step */
  double dt = (dx < dy) ? dx/c/50.0 : dy/c/50.0;

  /* start time stepping */
  double runtime = 0.0;
  double commtime = 0.0;
  for (it=0; it<nt; it++) {

    /* communicate to get v2 & v3 neighbor values */
    stime = MPI_Wtime();
    ierr = Communication1(v2, v3, &p2d);
    check_err(ierr, p2d.comm, "Communication1");
    ftime = MPI_Wtime();
    commtime += ftime - stime;

    /* start timer */
    stime = MPI_Wtime();

    /* first update v1 to get to half time step
       get relevant values for this location */
    for (j=0; j<p2d.nyloc; j++) {
      for (i=0; i<p2d.nxloc; i++) {

	/* access relevant components of v2 and v3 */
	v2_E = (i == p2d.nxloc-1) ? p2d.v2recvE[j] : v2[idx(i+1,j,p2d.nxloc)];
	v2_W = v2[idx(i,j,p2d.nxloc)];
	v3_N = (j == p2d.nyloc-1) ? p2d.v3recvN[i] : v3[idx(i,j+1,p2d.nxloc)];
	v3_S = v3[idx(i,j,p2d.nxloc)];

	/* update v1*/
	v1[idx(i,j,p2d.nxloc)] += c*dt/dx*(v2_E - v2_W)
                                + c*dt/dy*(v3_N - v3_S);

      } /* for i */
    } /* for j */
    ftime = MPI_Wtime();
    runtime += ftime - stime;

    /* communicate to get v1 neighbor values */
    stime = MPI_Wtime();
    ierr = Communication2(v1, &p2d);
    check_err(ierr, p2d.comm, "Communication2");
    ftime = MPI_Wtime();
    commtime += ftime - stime;

    /* next update v2 & v3 to get to full time step
       get relevant values for this location */
    stime = MPI_Wtime();
    for (j=0; j<p2d.nyloc; j++) {
      for (i=0; i<p2d.nxloc; i++) {

	/* access relevant components of v1 */
	v1_W = (i == 0) ? p2d.v1recvW[j] : v1[idx(i-1,j,p2d.nxloc)];
	v1_E = v1[idx(i,j,p2d.nxloc)];
	v1_S = (j == 0) ? p2d.v1recvS[i] : v1[idx(i,j-1,p2d.nxloc)];
	v1_N = v1[idx(i,j,p2d.nxloc)];

	/* update v2 and v3 */
	v2[idx(i,j,p2d.nxloc)] += c*dt/dx*(v1_E - v1_W);
	v3[idx(i,j,p2d.nxloc)] += c*dt/dy*(v1_N - v1_S);
	
      } /* for i */
    } /* for j */

    /* update solution for plotting */
    for (j=0; j<p2d.nyloc; j++) 
      for (i=0; i<p2d.nxloc; i++) 
	u[idx(i,j,p2d.nxloc)] += dt*v1[idx(i,j,p2d.nxloc)];

    /* update time */
    t = t + dt;
    ftime = MPI_Wtime();
    runtime += ftime - stime;

    /* stop simulation if we've reached tstop */
    if (t >= tstop)  break;

    /* output solution periodically */
    if ( t - (toutput+dtoutput) > -1.e-14) {
      stime = MPI_Wtime();
      toutput = t;
      noutput++;
      if (p2d.myid == 0)
	printf("writing output file %i, step = %i, t = %g\n", noutput, it, t);
      output(u, t, nx, ny, noutput, &p2d);
      ftime = MPI_Wtime();
      iotime += ftime - stime;
    }

  } /* for it */
  
  /* output final solution */
  stime = MPI_Wtime();
  toutput = t;
  noutput++;
  if (p2d.myid == 0)
    printf("writing output file %i, step = %i, t = %g\n", noutput, it, t);
  output(u, t, nx, ny, noutput, &p2d);
  ftime = MPI_Wtime();
  iotime += ftime - stime;

  if (p2d.myid == 0) {
    printf(" total initialization time = %.16e\n", inittime);
    printf(" total input/output time   = %.16e\n", iotime);
    printf(" total communication time  = %.16e\n", commtime);
    printf(" total simulation time     = %.16e\n", runtime);
  }

  /* free memory */
  free(u);
  free(v1);
  free(v2);
  free(v3);
  free_parallel_decomp(&p2d);

  /* finalize MPI */
  ierr = MPI_Finalize();

  return 0;
} /* end main */



void create_parallel_decomp(parallel_decomp *p2d) {
  p2d->pdims[0] = 1;
  p2d->pdims[1] = 1;
  p2d->periodic[0] = 1;
  p2d->periodic[1] = 1;
  p2d->pcoords[0] = 0;
  p2d->pcoords[1] = 0;
  p2d->nxloc = 0;
  p2d->nyloc = 0;
  p2d->myid = -1;
  p2d->numprocs = -1;
  p2d->nbE = -1;
  p2d->nbW = -1;
  p2d->nbN = -1;
  p2d->nbS = -1;
  p2d->v2recvE = NULL;
  p2d->v3recvN = NULL;
  p2d->v2sendW = NULL;
  p2d->v3sendS = NULL;
  p2d->v1recvW = NULL;
  p2d->v1recvS = NULL;
  p2d->v1sendE = NULL;
  p2d->v1sendN = NULL;
}



void free_parallel_decomp(parallel_decomp *p2d) {
  free(p2d->v2recvE);
  free(p2d->v3recvN);
  free(p2d->v2sendW);
  free(p2d->v3sendS);
  free(p2d->v1recvW);
  free(p2d->v1recvS);
  free(p2d->v1sendE);
  free(p2d->v1sendN);
}



int Communication1(double *v2, double *v3, parallel_decomp *p2d) {

  /* declarations */
  int i, j, ierr, nbcoords[2];
  
  /* allocate send/recv arrays if NULL */
  if (p2d->v2recvE == NULL)
    p2d->v2recvE = (double *) malloc( p2d->nyloc * sizeof(double) );
  if (p2d->v3recvN == NULL)
    p2d->v3recvN = (double *) malloc( p2d->nxloc * sizeof(double) );
  if (p2d->v2sendW == NULL)
    p2d->v2sendW = (double *) malloc( p2d->nyloc * sizeof(double) );
  if (p2d->v3sendS == NULL)
    p2d->v3sendS = (double *) malloc( p2d->nxloc * sizeof(double) );

  /* initialize send/recv buffers */
  for (j=0; j<p2d->nyloc; j++)  p2d->v2sendW[j] = v2[idx(0,j,p2d->nxloc)];
  for (i=0; i<p2d->nxloc; i++)  p2d->v3sendS[i] = v3[idx(i,0,p2d->nxloc)];
  for (j=0; j<p2d->nyloc; j++)  p2d->v2recvE[j] = 0.0;
  for (i=0; i<p2d->nxloc; i++)  p2d->v3recvN[i] = 0.0;

  /* determine process 'neighbors' */
  if (p2d->nbE < 0) {
    nbcoords[0] = p2d->pcoords[0]-1;
    nbcoords[1] = p2d->pcoords[1];
    ierr = MPI_Cart_rank(p2d->comm, nbcoords, &(p2d->nbW));
    check_err(ierr, p2d->comm, "MPI_Cart_rank");

    nbcoords[0] = p2d->pcoords[0]+1;
    nbcoords[1] = p2d->pcoords[1];
    ierr = MPI_Cart_rank(p2d->comm, nbcoords, &(p2d->nbE));
    check_err(ierr, p2d->comm, "MPI_Cart_rank");

    nbcoords[0] = p2d->pcoords[0];
    nbcoords[1] = p2d->pcoords[1]-1;
    ierr = MPI_Cart_rank(p2d->comm, nbcoords, &(p2d->nbS));
    check_err(ierr, p2d->comm, "MPI_Cart_rank");

    nbcoords[0] = p2d->pcoords[0];
    nbcoords[1] = p2d->pcoords[1]+1;
    ierr = MPI_Cart_rank(p2d->comm, nbcoords, &(p2d->nbN));
    check_err(ierr, p2d->comm, "MPI_Cart_rank");
  }  /* end if (nbE < 0) */

  /* open send/receive channels */
  MPI_Request req[4];
  ierr = MPI_Irecv(p2d->v2recvE, p2d->nyloc, MPI_DOUBLE,
                   p2d->nbE, 2, p2d->comm, &(req[0]));
  check_err(ierr, p2d->comm, "MPI_Irecv");

  ierr = MPI_Irecv(p2d->v3recvN, p2d->nxloc, MPI_DOUBLE,
                   p2d->nbN, 4, p2d->comm, &(req[1]));
  check_err(ierr, p2d->comm, "MPI_Irecv");

  ierr = MPI_Isend(p2d->v2sendW, p2d->nyloc, MPI_DOUBLE,
                   p2d->nbW, 2, p2d->comm, &(req[2]));
  check_err(ierr, p2d->comm, "MPI_Isend");

  ierr = MPI_Isend(p2d->v3sendS, p2d->nxloc, MPI_DOUBLE,
                   p2d->nbS, 4, p2d->comm, &(req[3]));
  check_err(ierr, p2d->comm, "MPI_Isend");

  /* wait until all communications finish */
  MPI_Status stat[4];
  ierr = MPI_Waitall(4, req, stat);
  check_err(ierr, p2d->comm, "MPI_Waitall");

  return 0;
}




int Communication2(double *v1, parallel_decomp *p2d) {
  
  /* declarations */
  int i, j, ierr, nbcoords[2];
  
  /* allocate send/recv arrays if NULL */
  if (p2d->v1recvW == NULL)
    p2d->v1recvW = (double *) malloc( p2d->nyloc * sizeof(double) );
  if (p2d->v1recvS == NULL)
    p2d->v1recvS = (double *) malloc( p2d->nxloc * sizeof(double) );
  if (p2d->v1sendE == NULL)
    p2d->v1sendE = (double *) malloc( p2d->nyloc * sizeof(double) );
  if (p2d->v1sendN == NULL)
    p2d->v1sendN = (double *) malloc( p2d->nxloc * sizeof(double) );

  /* initialize send/recv buffers */
  for (j=0; j<p2d->nyloc; j++)  p2d->v1sendE[j] = v1[idx(p2d->nxloc-1,j,p2d->nxloc)];
  for (i=0; i<p2d->nxloc; i++)  p2d->v1sendN[i] = v1[idx(i,p2d->nyloc-1,p2d->nxloc)];
  for (j=0; j<p2d->nyloc; j++)  p2d->v1recvW[j] = 0.0;
  for (i=0; i<p2d->nxloc; i++)  p2d->v1recvS[i] = 0.0;

  /* open send/receive channels */
  MPI_Request req[4];
  ierr = MPI_Irecv(p2d->v1recvW, p2d->nyloc, MPI_DOUBLE,
                   p2d->nbW, 2, p2d->comm, &(req[0]));
  check_err(ierr, p2d->comm, "MPI_Irecv");

  ierr = MPI_Irecv(p2d->v1recvS, p2d->nxloc, MPI_DOUBLE,
                   p2d->nbS, 4, p2d->comm, &(req[1]));
  check_err(ierr, p2d->comm, "MPI_Irecv");

  ierr = MPI_Isend(p2d->v1sendE, p2d->nyloc, MPI_DOUBLE,
                   p2d->nbE, 2, p2d->comm, &(req[2]));
  check_err(ierr, p2d->comm, "MPI_Isend");

  ierr = MPI_Isend(p2d->v1sendN, p2d->nxloc, MPI_DOUBLE,
                   p2d->nbN, 4, p2d->comm, &(req[3]));
  check_err(ierr, p2d->comm, "MPI_Isend");

  /* wait until all communications finish */
  MPI_Status stat[4];
  ierr = MPI_Waitall(4, req, stat);
  check_err(ierr, p2d->comm, "MPI_Waitall");

  return 0;
}
