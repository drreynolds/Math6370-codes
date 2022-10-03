cmake_minimum_required (VERSION 3.10)
project (Math6370ClassCodes)

set(Kokkos_DIR "$ENV{Kokkos_ROOT}" CACHE STRING "Kokkos root directory")
find_package(Kokkos REQUIRED)

add_executable(advection.exe advection.cpp initialize.cpp output.cpp)
target_link_libraries(advection.exe Kokkos::kokkos)

add_executable(advection_hybrid.exe advection_hybrid.cpp initialize.cpp output_hybrid.cpp)
target_link_libraries(advection.exe Kokkos::kokkos)
