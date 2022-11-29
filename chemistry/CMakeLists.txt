cmake_minimum_required (VERSION 3.10)
project (Math6370ClassCodes)

set(Kokkos_DIR "$ENV{Kokkos_ROOT}" CACHE STRING "Kokkos root directory")
find_package(Kokkos REQUIRED)

add_executable(chemistry.serial chemistry.cpp)
target_link_libraries(chemistry.serial Kokkos::kokkos)

add_executable(chemistry.openmp chemistry.cpp)
target_compile_definitions(chemistry.openmp PRIVATE -DUSE_OPENMP)
target_link_libraries(chemistry.openmp Kokkos::kokkos)

add_executable(chemistry.cuda chemistry.cpp)
target_compile_definitions(chemistry.cuda PRIVATE -DUSE_CUDA)
target_link_libraries(chemistry.cuda Kokkos::kokkos)

add_executable(chemistry.uvm chemistry.cpp)
target_compile_definitions(chemistry.uvm PRIVATE -DUSE_UVM)
target_link_libraries(chemistry.uvm Kokkos::kokkos)
