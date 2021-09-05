/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   7 February 2015 */

// Inclusions 
#include <stdlib.h>
#include <iostream>
#include "omp.h"

// OpenMP "schedule" example
int main(int argc, char* argv[]) {

  // declarations
  int i, n=23;

  // loop
  #pragma omp parallel for schedule(runtime)
  for (i=0; i<n; i++) 
    std::cout << "iteration" << i << " performed by thread " 
	      << omp_get_thread_num() << "\n";

  return 0;
} // end main 
