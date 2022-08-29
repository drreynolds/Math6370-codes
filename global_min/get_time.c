/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   1 February 2017 */


/* Inclusions */
#include "get_time.h"

/* Function to get the current "time" and cast to double.
   The raw returned value is unusable on its own; the utility 
   in this function comes when used as a 'stopwatch':
      double stime = get_time();
      ...
      double ftime = get_time();
      double elapsed_time = ftime-stime;
 */
double get_time() {
#ifdef _OPENMP
  return omp_get_wtime();
#else
  return ((1.0*clock()) / CLOCKS_PER_SEC);
#endif
}

/* Function to return the resolution of get_time */
double time_resolution() {
#ifdef _OPENMP
  return omp_get_wtick();
#else
  return (1.0 / CLOCKS_PER_SEC);
#endif
}
