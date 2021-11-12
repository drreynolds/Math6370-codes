/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   12 November 2021 */

// Inclusions
#include "advection.hpp"

// Writes current solution to disk
void output(double *u_h, double *u_d, double t, int nx, int ny, int noutput) {

  // copy device data to host
  cudaMemcpy( u_h, u_d, nx*ny*sizeof(double), cudaMemcpyDeviceToHost);

  // set output file name
  // Note: we reserve the first set of digits for the MPI process (unused here)
  char outname[100];
  sprintf(outname, "u_sol.000.%03i", noutput);

  // open output file
  FILE *FID = fopen(outname,"w");

  // write data set parameters
  // Note: the two 0's will be used for the MPI process location (unused here)
  fprintf(FID, "%i\n", nx);
  fprintf(FID, "%i\n", ny);
  fprintf(FID, "%i\n", 0);
  fprintf(FID, "%i\n", 0);
  fprintf(FID, "%.16e\n", t);

  // output the solution values and close the data set
  for (int j=0; j<ny; j++)
    for (int i=0; i<nx; i++)
      fprintf(FID, "%.16e\n",u_h[idx(i,j,nx)]);
  fclose(FID);

  // now output a metadata file, containing general run information
  // Note: the two 1's will be used for the MPI process dimensions (unused here)
  FID = fopen("u_sol.txt","w");
  fprintf(FID, "%i\n", nx);
  fprintf(FID, "%i\n", ny);
  fprintf(FID, "%i\n", 1);
  fprintf(FID, "%i\n", 1);
  fprintf(FID, "%i\n", noutput);
  fclose(FID);

} // end output
