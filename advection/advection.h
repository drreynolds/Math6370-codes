/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   7 February 2015 */

// simple macro to map a 2D index to a 1D address space
#define idx(i,j,nx)   ((j)*(nx)+(i))

// Prototypes
void initialize(double *u, double *v1, double *v2, double *v3, 
		double c, double dx, double dy, int nx, int ny);

void output(double *u, double t, int nx, int ny, int noutput);
