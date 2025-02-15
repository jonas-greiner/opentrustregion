cmake_minimum_required(VERSION 3.22)
project(trustorbopt Fortran)

# default build type
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")

# set flags for build types
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -Wextra -fcheck=all -fbacktrace")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -march=native")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -warn all -traceback")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -xHost")
endif()

# try finding 32-bit BLAS/LAPACK first
set(BLA_SIZEOF_INTEGER 4 CACHE STRING "Size of BLAS/LAPACK integer type (4 for LP64, 8 for ILP64)")
find_package(BLAS QUIET)
find_package(LAPACK QUIET)

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
add_library(trustorbopt SHARED src/trustorbopt.f90 src/c_interface.f90)

# enable preprocessing
set_source_files_properties(src/trustorbopt.f90 PROPERTIES Fortran_PREPROCESS ON)

# link against BLAS and LAPACK
target_link_libraries(trustorbopt PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# set position-independent code for shared libraries
set_target_properties(trustorbopt PROPERTIES POSITION_INDEPENDENT_CODE ON)
