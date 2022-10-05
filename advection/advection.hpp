/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// includes
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <Kokkos_Core.hpp>

// set the Kokkos execution space, memory spaces and range policies
#if defined(USE_OPENMP)
  typedef Kokkos::OpenMP     DevExec;
  typedef Kokkos::OpenMP     HostExec;
  typedef Kokkos::HostSpace  MemSpace;
#elif defined(USE_CUDA)
  typedef Kokkos::Cuda       DevExec;
  typedef Kokkos::Serial     HostExec;
  typedef Kokkos::CudaSpace  MemSpace;
#else
  typedef Kokkos::Serial     DevExec;
  typedef Kokkos::Serial     HostExec;
  typedef Kokkos::HostSpace  MemSpace;
#endif
typedef Kokkos::View<double**, MemSpace>  Vec2D;
typedef Vec2D::HostMirror                 Vec2DHost;
typedef Kokkos::MDRangePolicy<DevExec,  Kokkos::Rank<2>>  DevRange2D;
typedef Kokkos::MDRangePolicy<HostExec, Kokkos::Rank<2>>  HostRange2D;

// Prototypes
void initialize(Vec2DHost u_h, Vec2DHost v1_h, Vec2DHost v2_h, Vec2DHost v3_h,
                Vec2D u_d, Vec2D v1_d, Vec2D v2_d, Vec2D v3_d,
                double c, double dx, double dy, int nx, int ny);

void output(Vec2DHost u_h, Vec2D u_d, double t, int nx, int ny, int noutput);
