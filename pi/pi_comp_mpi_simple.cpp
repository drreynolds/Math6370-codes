/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "mpi.h"


// Prototypes
inline double f(double a) { return (4.0 / (1.0 + a*a)); }


/* Example routine to compute pi using numerical integration via
      pi = 4 * int_0^1 1/(1+x^2) dx
   We use a simple midpoint rule for integration, over
   subintervals of fixed size 1/n, where n is a user-input parameter. */
int main(int argc, char* argv[]) {

  // intialize MPI
  int ierr, numprocs, myid;
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // input the number of intervals
  int n;
  if (myid == 0) {
    std::cout << " Running with " << numprocs << " MPI tasks\n";
    std::cout << " Enter the number of intervals (0 quits):\n";
    std::cin >> n;
  }

  // root node sends n out to other processors
  ierr = MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // stop for illegal n
  if (n < 1)  MPI_Abort(MPI_COMM_WORLD, 1);

  // start timer
  double stime = MPI_Wtime();

  // set subinterval width
  const double h = 1.0 / n;

  // determine this processor's interval
  int is = ((int) (1.0*n/numprocs))*myid + 1;
  int ie = ((int) (1.0*n/numprocs))*(myid+1);
  if (myid == numprocs-1)  ie = n;

  // perform integration over n intervals
  double mypi = 0.0;
  for (int i=is; i<=ie; i++) {
    double x = h * (i - 0.5);
    mypi += h * f(x);
  }

  // root node collects result
  double pi;
  ierr = MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  // stop timer
  double ftime = MPI_Wtime();

  // output computed value and error
  if (myid == 0) {
    std::cout << " computed pi = " << std::setprecision(16) << pi << std::endl;
    std::cout << "     true pi = " << M_PI << std::endl;
    std::cout << "       error = " << M_PI-pi << std::endl;
    std::cout << "     runtime = " << ftime-stime << std::endl;
  }

  // finalize MPI
  ierr = MPI_Finalize();

}  // end main
