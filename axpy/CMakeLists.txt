cmake_minimum_required (VERSION 3.10)
project (Math6370ClassCodes)

set(Kokkos_DIR "$ENV{Kokkos_ROOT}" CACHE STRING "Kokkos root directory")
find_package(Kokkos REQUIRED)

add_executable(axpy.exe axpy.cpp)
target_link_libraries(axpy.exe Kokkos::kokkos)
