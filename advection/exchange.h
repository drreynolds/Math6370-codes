/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   7 February 2015 */

#ifndef _EXCHANGE_H
#define _EXCHANGE_H

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
  
  /* Set up an "exchange" structure, that will hold data buffers for 
     exchanges with neighboring processes.  Also prototype routines to:
        (a) initialize the structure
        (b) perform exchanges
        (c) free structure data. */

  // exchange_data structure definition
  typedef struct {
    int nxloc;         // local subdomain size
    int nyloc;
    int nbN;           // process IDs for nearest neighbors
    int nbS;
    int nbE;
    int nbW;
    MPI_Comm comm;     // communicator
    double *send_EW;   // send buffers
    double *send_NS;
  } exchange_data;

  // Function prototypes for acting on exchange structure

  //   creates the exchange_data structure
  void init_exchange(MPI_Comm com, int nxl, int nyl, 
		     int pcoords[2], exchange_data *exchange);

  //   frees memory from the exchange_data structure
  void free_exchange(exchange_data *exchange);

  //   performs exchange on v1 field
  int do_exchange1(double *v1, double *v1W, double *v1S, 
		   exchange_data *exchange);

  //   performs exchange on v2 and v3 fields
  int do_exchange2(double *v2, double *v3, double *v2E, 
		   double *v3N, exchange_data *exchange);

#endif


/*********************** END OF FILE ***********************/
