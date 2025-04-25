cmake_minimum_required (VERSION 3.10)
project (Math6370ClassCodes)

set(Kokkos_DIR "$ENV{Kokkos_ROOT}" CACHE STRING "Kokkos root directory")
find_package(Kokkos REQUIRED)

add_executable(advection.serial advection.cpp initialize.cpp output.cpp)
target_link_libraries(advection.serial Kokkos::kokkos)

add_executable(advection.cuda advection.cpp initialize.cpp output.cpp)
target_compile_definitions(advection.cuda PRIVATE -DUSE_CUDA)
target_link_libraries(advection.cuda Kokkos::kokkos)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
   add_executable(advection.hybrid advection_hybrid.cpp initialize.cpp output_hybrid.cpp)
   target_compile_definitions(advection.hybrid PRIVATE -DUSE_CUDA)
   target_link_libraries(advection.hybrid PUBLIC OpenMP::OpenMP_CXX Kokkos::kokkos)
endif()

