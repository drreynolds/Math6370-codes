cmake_minimum_required (VERSION 3.10)
project (Math6370ClassCodes)

set(Kokkos_DIR "$ENV{Kokkos_ROOT}" CACHE STRING "Kokkos root directory")
find_package(Kokkos REQUIRED)

add_executable(glob_min.serial glob_min.cpp)
target_link_libraries(glob_min.serial Kokkos::kokkos)

add_executable(glob_min.openmp glob_min.cpp)
target_compile_definitions(glob_min.openmp PRIVATE -DUSE_OPENMP)
target_link_libraries(glob_min.openmp Kokkos::kokkos)

add_executable(glob_min.cuda glob_min.cpp)
target_compile_definitions(glob_min.cuda PRIVATE -DUSE_CUDA)
target_link_libraries(glob_min.cuda Kokkos::kokkos)
