/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */


// Inclusions
#include <stdio.h>
#include <string.h>
#include "advection_mpi.hpp"


// Writes current solution to disk
void output(double* u, double t, int nx, int ny, int noutput, parallel_decomp& p2d) {

  // set process-local output file name
  char outname[100];
  sprintf(outname, "u_sol.%03i.%03i", p2d.myid, noutput);

  // open process-local output file
  FILE *FID = fopen(outname,"w");

  // write data set parameters
  fprintf(FID, "%i\n", p2d.nxloc);
  fprintf(FID, "%i\n", p2d.nyloc);
  fprintf(FID, "%i\n", p2d.pcoords[0]);
  fprintf(FID, "%i\n", p2d.pcoords[1]);
  fprintf(FID, "%.16e\n", t);

  // output the solution values and close the data set
  for (int j=0; j<p2d.nyloc; j++)
    for (int i=0; i<p2d.nxloc; i++)
      fprintf(FID, "%.16e\n",u[idx(i,j,p2d.nxloc)]);
  fclose(FID);

  // root process outputs a metadata file, containing general run information
  if (p2d.myid == 0) {
    FID = fopen("u_sol.txt","w");
    fprintf(FID, "%i\n", nx);
    fprintf(FID, "%i\n", ny);
    fprintf(FID, "%i\n", p2d.pdims[0]);
    fprintf(FID, "%i\n", p2d.pdims[1]);
    fprintf(FID, "%i\n", noutput);
    fclose(FID);
  }

} // end output
