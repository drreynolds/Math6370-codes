cmake_minimum_required (VERSION 3.10)
project (Math6370ClassCodes)

set(Kokkos_DIR "$ENV{Kokkos_ROOT}" CACHE STRING "Kokkos root directory")
find_package(Kokkos REQUIRED)

add_executable(dot_prod.serial dot_prod.cpp)
target_link_libraries(dot_prod.serial Kokkos::kokkos)

add_executable(dot_prod.openmp dot_prod.cpp)
target_compile_definitions(dot_prod.openmp PRIVATE -DUSE_OPENMP)
target_link_libraries(dot_prod.openmp Kokkos::kokkos)

add_executable(dot_prod.cuda dot_prod.cpp)
target_compile_definitions(dot_prod.cuda PRIVATE -DUSE_CUDA)
target_link_libraries(dot_prod.cuda Kokkos::kokkos)
