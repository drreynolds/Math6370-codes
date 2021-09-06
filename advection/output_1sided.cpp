/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   30 April 2015 */


// Inclusions
#include <stdio.h>
#include <string.h>
#include "advection_1sided.h"


// Writes current solution to disk
void output(double* u, double t, int nx, int ny, int nxl, 
	    int nyl, int noutput, MPI_Comm comm) {

  // declarations
  char outname[100];
  FILE* FID;
  int i, j, myid, ierr, myloc[2], pdims[2], periods[2];

  // get myid from comm
  ierr = MPI_Comm_rank(comm, &myid);
  check_err(ierr, comm, "MPI_Comm_rank");

  // get process location from comm
  ierr = MPI_Cart_get(comm, 2, pdims, periods, myloc);
  check_err(ierr, comm, "MPI_Cart_get");

  // set output file name
  sprintf(outname, "u_sol.%03i.%03i", myid, noutput);

  // open output file
  FID = fopen(outname,"w");

  // write data set parameters
  fprintf(FID, "%i\n", nxl);
  fprintf(FID, "%i\n", nyl);
  fprintf(FID, "%i\n", myloc[0]);
  fprintf(FID, "%i\n", myloc[1]);
  fprintf(FID, "%.16e\n", t);

  // output the solution values and close the data set
  for (j=0; j<nyl; j++) 
    for (i=0; i<nxl; i++) 
      fprintf(FID, "%.16e\n",u[idx(i,j,nxl)]);
  fclose(FID);
    
  // now output a metadata file, containing general run information
  if (myid == 0) {
    FID = fopen("u_sol.txt","w");
    fprintf(FID, "%i\n", nx);
    fprintf(FID, "%i\n", ny);
    fprintf(FID, "%i\n", pdims[0]);
    fprintf(FID, "%i\n", pdims[1]);
    fprintf(FID, "%i\n", noutput);
    fclose(FID);
  }

} // end output
