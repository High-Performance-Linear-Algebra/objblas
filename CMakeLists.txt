cmake_minimum_required(VERSION 3.28)

project(objblas LANGUAGES Fortran)

find_package(BLAS REQUIRED)
find_package(OpenMP REQUIRED COMPONENTS Fortran)

set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -mtune=native")
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)

add_library(objblas src/blas.f90 src/blasOMP.f90)
target_link_libraries(objblas PUBLIC BLAS::BLAS OpenMP::OpenMP_Fortran)

add_executable(test_axpy test/test_axpy.f90)
target_link_libraries(test_axpy PRIVATE objblas)

add_executable(axpy_perf test/axpy_perf.f90)
target_link_libraries(axpy_perf PRIVATE objblas)

add_executable(axpy_scaling test/axpy_scaling.f90)
target_link_libraries(axpy_scaling PRIVATE objblas)