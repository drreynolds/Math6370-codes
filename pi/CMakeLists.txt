cmake_minimum_required (VERSION 3.10)
project (Math6370ClassCodes)

set(Kokkos_DIR "$ENV{Kokkos_ROOT}" CACHE STRING "Kokkos root directory")
find_package(Kokkos REQUIRED)

add_executable(pi_comp.serial pi_comp.cpp)
target_link_libraries(pi_comp.serial Kokkos::kokkos)

add_executable(pi_comp.openmp pi_comp.cpp)
target_compile_definitions(pi_comp.openmp PRIVATE -DUSE_OPENMP)
target_link_libraries(pi_comp.openmp Kokkos::kokkos)

add_executable(pi_comp.cuda pi_comp.cpp)
target_compile_definitions(pi_comp.cuda PRIVATE -DUSE_CUDA)
target_link_libraries(pi_comp.cuda Kokkos::kokkos)
