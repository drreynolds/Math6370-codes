/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   7 February 2015 */

/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include "omp.h"

/* OpenMP "schedule" example */
int main(int argc, char* argv[]) {

  /* declarations */
  int i, n=23;

  /* loop */
  #pragma omp parallel for schedule(runtime)
  for (i=0; i<n; i++) 
    printf("iteration %i performed by thread %i\n", i, omp_get_thread_num());

  return 0;
} /* end main */
