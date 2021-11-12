/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   7 February 2015 */

// includes
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <RAJA/RAJA.hpp>
#include <string.h>

// simple macro to map a 2D index to a 1D address space
#define idx(i,j,nx)   ((j)*(nx)+(i))

// Prototypes
void initialize(double *u_h, double *v1_h, double *v2_h, double *v3_h, double *u_d, double *v1_d, double *v2_d, double *v3_d, double c, double dx, double dy, int nx, int ny);

void output(double *u_h, double *u_d, double t, int nx, int ny, int noutput);
