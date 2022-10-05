cmake_minimum_required (VERSION 3.10)
project (Math6370ClassCodes)

set(Kokkos_DIR "$ENV{Kokkos_ROOT}" CACHE STRING "Kokkos root directory")
find_package(Kokkos REQUIRED)

add_executable(advection.serial advection.cpp initialize.cpp output.cpp)
target_link_libraries(advection.serial Kokkos::kokkos)

add_executable(advection.openmp advection.cpp initialize.cpp output.cpp)
target_compile_definitions(advection.openmp PRIVATE -DUSE_OPENMP)
target_link_libraries(advection.openmp Kokkos::kokkos)

add_executable(advection.cuda advection.cpp initialize.cpp output.cpp)
target_compile_definitions(advection.cuda PRIVATE -DUSE_CUDA)
target_link_libraries(advection.cuda Kokkos::kokkos)

# add_executable(advection_hybrid.serial advection_hybrid.cpp initialize.cpp output_hybrid.cpp)
# target_link_libraries(advection_hybrid.serial Kokkos::kokkos)

# add_executable(advection_hybrid.openmp advection_hybrid.cpp initialize.cpp output_hybrid.cpp)
# target_compile_definitions(advection_hybrid.openmp PRIVATE -DUSE_OPENMP)
# target_link_libraries(advection_hybrid.openmp Kokkos::kokkos)

# add_executable(advection_hybrid.cuda advection_hybrid.cpp initialize.cpp output_hybrid.cpp)
# target_compile_definitions(advection_hybrid.cuda PRIVATE -DUSE_CUDA)
# target_link_libraries(advection_hybrid.cuda Kokkos::kokkos)
