cmake_minimum_required (VERSION 3.10)
project (Math6370ClassCodes)

set(Kokkos_DIR "$ENV{Kokkos_ROOT}" CACHE STRING "Kokkos root directory")
find_package(Kokkos REQUIRED)

add_executable(axpy.serial axpy.cpp)
target_link_libraries(axpy.serial Kokkos::kokkos)

add_executable(axpy.openmp axpy.cpp)
target_compile_definitions(axpy.openmp PRIVATE -DUSE_OPENMP)
target_link_libraries(axpy.openmp Kokkos::kokkos)

add_executable(axpy.cuda axpy.cpp)
target_compile_definitions(axpy.cuda PRIVATE -DUSE_CUDA)
target_link_libraries(axpy.cuda Kokkos::kokkos)
