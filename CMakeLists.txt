cmake_minimum_required(VERSION 3.22)
project(opentrustregion Fortran C)

# default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# set flags for build types
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -Wextra -fcheck=all -fbacktrace")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -march=native")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "-fpscomp logicals")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -warn all -traceback")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -xHost")
endif()

# try finding 32-bit BLAS/LAPACK first
set(BLA_SIZEOF_INTEGER 4 CACHE STRING "Size of BLAS/LAPACK integer type (4 for LP64, 8 for ILP64)")
find_package(BLAS)
find_package(LAPACK)

# check if 32-bit BLAS/LAPACK found
if (NOT BLAS_FOUND OR NOT LAPACK_FOUND)
    message(STATUS "BLAS/LAPACK with 32-bit integers not found, trying version with 64-bit integers...")

    # now try finding 64-bit BLAS/LAPACK
    set(BLA_SIZEOF_INTEGER 8 CACHE STRING "Size of BLAS/LAPACK integer type (4 for LP64, 8 for ILP64)" FORCE)
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)

    if (BLAS_FOUND AND LAPACK_FOUND)
        message(STATUS "Using BLAS/LAPACK with 64-bit integers (ILP64) ")
        add_definitions(-DUSE_ILP64)
    else()
        message(FATAL_ERROR "Could not find a suitable BLAS/LAPACK library (neither with 32-bit or 64-bit integers).")
    endif()
else()
    message(STATUS "Using BLAS/LAPACK with 32-bit integers (LP64)")
endif()

# add source files
add_library(opentrustregion SHARED src/opentrustregion.f90 src/c_interface.f90)

# enable preprocessing
set_source_files_properties(src/opentrustregion.f90 PROPERTIES Fortran_PREPROCESS ON)

# link against BLAS and LAPACK
target_link_libraries(opentrustregion PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# create test executable
add_library(
    testsuite SHARED
    tests/opentrustregion_unit_tests.f90
    tests/c_interface_unit_tests.f90
    tests/opentrustregion_mock.f90
    tests/c_interface_mock.f90
    tests/opentrustregion_system_tests.f90
)
target_link_libraries(
    testsuite 
    PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} opentrustregion
)

# get the absolute path to the test data directory and add it as a definition so the
# Fortran code can access it
set(TEST_DATA_PATH "${CMAKE_SOURCE_DIR}/tests/data")
add_definitions(-DTEST_DATA_PATH=\"${TEST_DATA_PATH}\")
target_compile_definitions(testsuite PRIVATE TEST_DATA_PATH="${TEST_DATA_PATH}")

# enable preprocessing
set_source_files_properties(tests/opentrustregion_system_tests.f90 PROPERTIES Fortran_PREPROCESS ON)

# set paths to find module files and library in build directory
set(OpenTrustRegion_INCLUDE_DIR ${CMAKE_BINARY_DIR})
set(OpenTrustRegion_LIBRARY ${CMAKE_BINARY_DIR}/libopentrustregion${CMAKE_SHARED_LIBRARY_SUFFIX})

# create the OpenTrustRegionConfig.cmake from the template
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/OpenTrustRegionConfig.cmake.in
  ${CMAKE_BINARY_DIR}/cmake/OpenTrustRegionConfig.cmake
  @ONLY
)
