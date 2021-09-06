/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   30 April 2015 */


// Inclusions
#include <stdio.h>
#include <string.h>
#include "mpi.h"


// error checking routine for successful MPI calls
void check_err(const int ierr, MPI_Comm comm, const char* fname) {
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr, " error in %s = %i\n", fname, ierr);
    MPI_Abort(comm, ierr);
  }
}
