cmake_minimum_required(VERSION 3.24)
project(OpenTrustRegion VERSION 3.0.0 LANGUAGES Fortran C)

# get commands
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CheckFortranFunctionExists)

# set namespace
set(otr OpenTrustRegion)

# set components
set(runtime_component runtime)
set(test_component testsuite)
set(python_component python)

# ====  Options  ================================================================

# set installation directories, this needs to be separate for windows and unix to avoid 
# repeating ${otr} twice in the path
if ("${CMAKE_SYSTEM_NAME}" STREQUAL "Windows")
    set(${otr}_INSTALL_CMAKEDIR "${otr}/CMake"
        CACHE STRING "Directory to which CMake files are installed")
else()
    set(${otr}_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/${otr}"
        CACHE STRING "Directory to which CMake files are installed")
endif()
message(STATUS "Showing option ${otr}_INSTALL_CMAKEDIR: ${${otr}_INSTALL_CMAKEDIR}")

# build the testsuite by default only as the top-level project
option(${otr}_BUILD_TESTING "Build testsuite" ${PROJECT_IS_TOP_LEVEL})
option(${otr}_INSTALL_TESTING_LIBRARY "Install the compiled library used for testsuite" OFF)
message(STATUS "Showing option ${otr}_BUILD_TESTING: ${${otr}_BUILD_TESTING}")
message(STATUS "Showing option ${otr}_INSTALL_TESTING_LIBRARY: ${${otr}_INSTALL_TESTING_LIBRARY}")

# host program provides BLAS/LAPACK (skip find/link in this project)
option(${otr}_HOST_PROVIDES_BLAS "Host program will provide BLAS/LAPACK; skip finding/linking them" OFF)
message(STATUS "Showing option ${otr}_HOST_PROVIDES_BLAS: ${${otr}_HOST_PROVIDES_BLAS}")

# default to static libraries
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
message(STATUS "Showing option BUILD_SHARED_LIBS: ${BUILD_SHARED_LIBS}")

option(${otr}_ENABLE_XHOST "Optimize for current machine's instruction set" OFF)
message(STATUS "Showing option ${otr}_ENABLE_XHOST: ${${otr}_ENABLE_XHOST}")

# optional extension modules, all off by default
option(ENABLE_MO "Include molecular orbital extension in OpenTrustRegion" OFF)
option(ENABLE_OAO "Include orthogonal AO extension in OpenTrustRegion" OFF)
option(ENABLE_ROTATION_QN "Include rotation quasi-Newton extension in OpenTrustRegion" OFF)
option(ENABLE_DENSITY_QN "Include density quasi-Newton extension in OpenTrustRegion" OFF)
option(ENABLE_S_GEK
       "Include subspace gradient-enhanced kriging extension in OpenTrustRegion" OFF)
# density quasi-Newton works on the objects of both orbital bases, which it switches on
# for this configuration only, so that they are off again once it is switched off
if(ENABLE_DENSITY_QN AND NOT ENABLE_MO)
    message(STATUS "Including molecular orbital extension module as a dependency for density quasi-Newton")
    set(ENABLE_MO ON)
endif()
if(ENABLE_DENSITY_QN AND NOT ENABLE_OAO)
    message(STATUS "Including orthogonal AO extension module as a dependency for density quasi-Newton")
    set(ENABLE_OAO ON)
endif()
message(STATUS "Showing option ENABLE_MO: ${ENABLE_MO}")
message(STATUS "Showing option ENABLE_OAO: ${ENABLE_OAO}")
message(STATUS "Showing option ENABLE_ROTATION_QN: ${ENABLE_ROTATION_QN}")
message(STATUS "Showing option ENABLE_DENSITY_QN: ${ENABLE_DENSITY_QN}")
message(STATUS "Showing option ENABLE_S_GEK: ${ENABLE_S_GEK}")

# every extension builds on the code the extensions share
if(ENABLE_MO OR ENABLE_OAO OR ENABLE_ROTATION_QN OR ENABLE_DENSITY_QN OR ENABLE_S_GEK)
    set(ANY_EXTENSION_ENABLED ON)
else()
    set(ANY_EXTENSION_ENABLED OFF)
endif()

# ====  Build  ==================================================================

# default build type of a top-level build, a parent project that adds this one as a
# subdirectory chooses its own
if(PROJECT_IS_TOP_LEVEL AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# set flags for build types
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -Wextra -fcheck=all -fbacktrace")

    # compile every procedure for recursion in every build type, so that a call can be
    # nested in a callback function of another call or run concurrently in a thread
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-frecursive>")

    if(CMAKE_CROSSCOMPILING OR NOT ${otr}_ENABLE_XHOST)
        set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
    else()
        set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -march=native")
    endif()

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -warn all -traceback")

    # store .true. as 1 rather than -1 in every build type, since 1 is the only true
    # value of the C bool that logical(c_bool) is passed to
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-fpscomp logicals>")

    # compile every procedure for recursion in every build type, so that a call can be
    # nested in a callback function of another call or run concurrently in a thread
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-recursive>")

    if(CMAKE_CROSSCOMPILING OR NOT ${otr}_ENABLE_XHOST)
        set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
    else()
        set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -xHost")
    endif()

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "LLVMFlang")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g")

    # Flang compiles every procedure for recursion by default, as Fortran 2018 asks, so
    # nested and concurrent calls need no flag

    if(CMAKE_CROSSCOMPILING OR NOT ${otr}_ENABLE_XHOST)
        set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
    else()
        set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -march=native")
    endif()

    # CMake has no Apple platform file for Flang and passes the C compiler's
    # -dynamiclib (before CMake 4.0) and bare -install_name, which only the C compiler
    # driver accepts
    if(APPLE)
        set(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS "-shared")
        set(CMAKE_SHARED_LIBRARY_SONAME_Fortran_FLAG "-Xlinker -install_name -Xlinker")
    endif()

elseif(NOT CMAKE_Fortran_COMPILER_ID)
    # CMake sets no flags for a compiler it does not recognize, so the library would
    # be built without optimization
    message(WARNING "Fortran compiler ${CMAKE_Fortran_COMPILER} not recognized by "
                    "CMake ${CMAKE_VERSION}, so no optimization flags are set; pass "
                    "them through CMAKE_Fortran_FLAGS or use a newer CMake")
endif()

# user-supplied integer size
set(INTEGER_SIZE "" CACHE STRING "Size of integers for library (4=32-bit, 8=64-bit)")
if(INTEGER_SIZE AND NOT INTEGER_SIZE MATCHES "^(4|8)$")
    message(FATAL_ERROR
        "INTEGER_SIZE is ${INTEGER_SIZE}. "
        "Please set INTEGER_SIZE to 4 (LP64) or 8 (ILP64), or leave it unset to autodetect."
    )
endif()

# user-supplied BLAS and LAPACK libraries
set(BLAS_LIBRARIES "" CACHE STRING "User-specified BLAS library path(s)")
set(LAPACK_LIBRARIES "" CACHE STRING "User-specified LAPACK library path(s)")

# check if integer size is provided
if(INTEGER_SIZE)
    set(BLA_SIZEOF_INTEGER ${INTEGER_SIZE} CACHE INTERNAL "BLAS/LAPACK integer size")
    if(NOT ${otr}_HOST_PROVIDES_BLAS)
        if(BLAS_LIBRARIES)
            message(STATUS "Using user-specified BLAS")
        else()
            find_package(BLAS REQUIRED)
        endif()
        if(LAPACK_LIBRARIES)
            message(STATUS "Using user-specified LAPACK")
        else()
            find_package(LAPACK REQUIRED)
        endif()
    else()
        message(STATUS 
            "Host will provide BLAS/LAPACK; skipping finding and linking these libraries "
            "and building testsuite since tests require BLAS at build/link time"
        )
    endif()
else()
    if(BLAS_LIBRARIES OR LAPACK_LIBRARIES)
        message(FATAL_ERROR
            "Explicitly provided BLAS or LAPACK libraries require INTEGER_SIZE to be set. "
            "Please set INTEGER_SIZE to 4 (LP64) or 8 (ILP64) to ensure compatibility."
        )
    elseif(${otr}_HOST_PROVIDES_BLAS)
        message(FATAL_ERROR
            "OpenTrustRegion_HOST_PROVIDES_BLAS requires INTEGER_SIZE to be set. "
            "Please set INTEGER_SIZE to 4 (LP64) or 8 (ILP64) to match the host's BLAS/LAPACK integer width."
        )
    else()
        # autodetect BLAS and LAPACK, try 32-bit first, then 64-bit
        set(BLA_SIZEOF_INTEGER 4 CACHE INTERNAL "BLAS/LAPACK integer size")
        find_package(BLAS)
        find_package(LAPACK)
        if(NOT BLAS_FOUND OR NOT LAPACK_FOUND)
            set(BLA_SIZEOF_INTEGER 8 CACHE INTERNAL "BLAS/LAPACK integer size")
            find_package(BLAS REQUIRED)
            find_package(LAPACK REQUIRED)
        endif()
    endif()
endif()

# set integer size and define library suffix
message(STATUS "Using integer size ${BLA_SIZEOF_INTEGER}")
if(BLA_SIZEOF_INTEGER EQUAL 8)
    set(LIB_SUFFIX "64")

    # check whether 64-bit BLAS and LAPACK routines are called differently
    if(NOT ${otr}_HOST_PROVIDES_BLAS)
        set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
        check_fortran_function_exists("sgemm_64" BLAS_64)
        check_fortran_function_exists("sgesv_64" LAPACK_64)
        if(BLAS_64)
            add_compile_definitions("ddot=ddot_64")
            add_compile_definitions("dnrm2=dnrm2_64")
            add_compile_definitions("dgemv=dgemv_64")
            add_compile_definitions("dgemm=dgemm_64")
            add_compile_definitions("dtrsm=dtrsm_64")
            add_compile_definitions("dsyrk=dsyrk_64")
        endif()
        if(LAPACK_64)
            add_compile_definitions("dsyev=dsyev_64")
            add_compile_definitions("dgeev=dgeev_64")
            add_compile_definitions("dgelss=dgelss_64")
            add_compile_definitions("dgesv=dgesv_64")
            add_compile_definitions("dgetrf=dgetrf_64")
            add_compile_definitions("dgetrs=dgetrs_64")
            add_compile_definitions("dposv=dposv_64")
            add_compile_definitions("dpstrf=dpstrf_64")
            add_compile_definitions("dpttrf=dpttrf_64")
            add_compile_definitions("zheev=zheev_64")
        endif()
    else()
        set(BLAS_64 FALSE)
        set(LAPACK_64 FALSE)
    endif()
else()
    set(LIB_SUFFIX "32")
endif()

# define base source files
set(OPENTRUSTREGION_SOURCES src/opentrustregion.f90 src/c_interface.f90)

# base include directories
set(OPENTRUSTREGION_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)

# conditionally add extension files, the code the extensions share once
if(ANY_EXTENSION_ENABLED)
    list(APPEND OPENTRUSTREGION_SOURCES
        extensions/common/src/common.f90
        extensions/common/src/common_c_interface.f90
    )
    list(APPEND OPENTRUSTREGION_INCLUDE_DIRS
        ${CMAKE_CURRENT_SOURCE_DIR}/extensions/common/include
    )
endif()
if(ENABLE_MO)
    message(STATUS "Including molecular orbitals in OpenTrustRegion")
    list(APPEND OPENTRUSTREGION_SOURCES
        extensions/mo/src/mo.f90
        extensions/mo/src/mo_c_interface.f90
    )
    list(APPEND OPENTRUSTREGION_INCLUDE_DIRS
        ${CMAKE_CURRENT_SOURCE_DIR}/extensions/mo/include
    )
endif()
if(ENABLE_OAO)
    message(STATUS "Including orthogonal AO in OpenTrustRegion")
    list(APPEND OPENTRUSTREGION_SOURCES
        extensions/oao/src/oao.f90
        extensions/oao/src/oao_c_interface.f90
    )
    list(APPEND OPENTRUSTREGION_INCLUDE_DIRS
        ${CMAKE_CURRENT_SOURCE_DIR}/extensions/oao/include
    )
endif()
if(ENABLE_ROTATION_QN)
    message(STATUS "Including rotation quasi-Newton in OpenTrustRegion")
    list(APPEND OPENTRUSTREGION_SOURCES
        extensions/rotation_qn/src/rotation_qn.f90
        extensions/rotation_qn/src/rotation_qn_c_interface.f90
    )
    list(APPEND OPENTRUSTREGION_INCLUDE_DIRS
        ${CMAKE_CURRENT_SOURCE_DIR}/extensions/rotation_qn/include
    )
endif()
if(ENABLE_DENSITY_QN)
    message(STATUS "Including density quasi-Newton in OpenTrustRegion")
    list(APPEND OPENTRUSTREGION_SOURCES
        extensions/density_qn/src/density_qn.f90
        extensions/density_qn/src/density_qn_c_interface.f90
    )
    list(APPEND OPENTRUSTREGION_INCLUDE_DIRS
        ${CMAKE_CURRENT_SOURCE_DIR}/extensions/density_qn/include
    )
endif()
if(ENABLE_S_GEK)
    message(STATUS "Including subspace gradient-enhanced kriging in OpenTrustRegion")
    list(APPEND OPENTRUSTREGION_SOURCES
        extensions/s_gek/src/s_gek.f90
        extensions/s_gek/src/s_gek_c_interface.f90
    )
    list(APPEND OPENTRUSTREGION_INCLUDE_DIRS
        ${CMAKE_CURRENT_SOURCE_DIR}/extensions/s_gek/include
    )
endif()

# add library
add_library(opentrustregion ${OPENTRUSTREGION_SOURCES})
add_library(${otr}::opentrustregion ALIAS opentrustregion)

# the library's Fortran sources select their integer kinds with USE_ILP64, while C
# programs compiled against the header of a 64-bit integer library need the same
# integer width through the namespaced OTR_ILP64, which is the only definition
# exported with the library target
if(BLA_SIZEOF_INTEGER EQUAL 8)
    target_compile_definitions(opentrustregion PRIVATE USE_ILP64 INTERFACE OTR_ILP64)
endif()

# set properties
set_target_properties(
    opentrustregion PROPERTIES
    # output name
    OUTPUT_NAME "opentrustregion_${LIB_SUFFIX}"
    EXPORT_NAME "opentrustregion${LIB_SUFFIX}"
    # library version
    VERSION ${PROJECT_VERSION}
    # interface version, which includes the minor version since a minor version can
    # change the settings and state types, so that a program built against another
    # minor version does not load the library
    SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    # position-independent code - not automatically enabled for static libraries, but
    # required when linking a static library into a shared library
    POSITION_INDEPENDENT_CODE ON
)

# allow undefined symbols at link time on macOS only when the host program supplies
# BLAS/LAPACK at runtime
if(${BUILD_SHARED_LIBS} AND APPLE AND ${otr}_HOST_PROVIDES_BLAS)
    set_target_properties(opentrustregion PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup")
endif()

# create target-specific version property
set(export_properties
    "${otr}_VERSION"
)
set_target_properties(opentrustregion PROPERTIES ${otr}_VERSION ${${otr}_VERSION})
set_property(TARGET opentrustregion APPEND PROPERTY EXPORT_PROPERTIES "${export_properties}")

# enable preprocessing
set_source_files_properties(${OPENTRUSTREGION_SOURCES}
                            PROPERTIES Fortran_PREPROCESS ON
)

# link against BLAS and LAPACK (unless host provides these)
if(NOT ${otr}_HOST_PROVIDES_BLAS)
    target_link_libraries(opentrustregion
                          PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
    )
endif()

# set include directories
set(module_dir ${CMAKE_CURRENT_BINARY_DIR}/modules)
set(module_install_dir ${CMAKE_INSTALL_INCLUDEDIR}/opentrustregion_${LIB_SUFFIX})
set_target_properties(opentrustregion PROPERTIES Fortran_MODULE_DIRECTORY ${module_dir})
target_include_directories(opentrustregion PUBLIC
                           "$<BUILD_INTERFACE:${OPENTRUSTREGION_INCLUDE_DIRS}>"
                           $<BUILD_INTERFACE:${module_dir}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
                           $<INSTALL_INTERFACE:${module_install_dir}>
)

# option to build testsuite but skip when the host provides BLAS/LAPACK (tests require 
# BLAS at build/link time)
if (${otr}_BUILD_TESTING AND NOT ${otr}_HOST_PROVIDES_BLAS)
    # define test files
    set(OPENTRUSTREGION_TESTS
        tests/test_reference.f90
        tests/opentrustregion_unit_tests.f90
        tests/c_interface_unit_tests.f90
        tests/opentrustregion_mock.f90
        tests/c_interface_mock.f90
        tests/opentrustregion_system_tests.f90
        tests/c_integration_tests.c)

    # conditionally add extension tests, those of the code the extensions share once
    if(ANY_EXTENSION_ENABLED)
        list(APPEND OPENTRUSTREGION_TESTS
            extensions/common/tests/common_mock.f90
            extensions/common/tests/common_test_reference.f90
            extensions/common/tests/common_unit_tests.f90
            extensions/common/tests/common_c_interface_unit_tests.f90
            extensions/common/tests/common_c_interface_mock.f90
            extensions/common/tests/common_c_integration_tests.c
        )
    endif()
    if(ENABLE_MO)
        list(APPEND OPENTRUSTREGION_TESTS
            extensions/mo/tests/mo_test_reference.f90
            extensions/mo/tests/mo_unit_tests.f90
            extensions/mo/tests/mo_c_interface_unit_tests.f90
            extensions/mo/tests/mo_mock.f90
            extensions/mo/tests/mo_c_interface_mock.f90
            extensions/mo/tests/mo_c_integration_tests.c
        )
    endif()
    if(ENABLE_OAO)
        list(APPEND OPENTRUSTREGION_TESTS
            extensions/oao/tests/oao_test_reference.f90
            extensions/oao/tests/oao_unit_tests.f90
            extensions/oao/tests/oao_c_interface_unit_tests.f90
            extensions/oao/tests/oao_mock.f90
            extensions/oao/tests/oao_c_interface_mock.f90
            extensions/oao/tests/oao_c_integration_tests.c
        )
    endif()
    if(ENABLE_ROTATION_QN)
        list(APPEND OPENTRUSTREGION_TESTS
            extensions/rotation_qn/tests/rotation_qn_test_reference.f90
            extensions/rotation_qn/tests/rotation_qn_unit_tests.f90
            extensions/rotation_qn/tests/rotation_qn_c_interface_unit_tests.f90
            extensions/rotation_qn/tests/rotation_qn_mock.f90
            extensions/rotation_qn/tests/rotation_qn_c_interface_mock.f90
            extensions/rotation_qn/tests/rotation_qn_c_integration_tests.c
        )
    endif()
    if(ENABLE_DENSITY_QN)
        list(APPEND OPENTRUSTREGION_TESTS
            extensions/density_qn/tests/density_qn_test_reference.f90
            extensions/density_qn/tests/density_qn_unit_tests.f90
            extensions/density_qn/tests/density_qn_c_interface_unit_tests.f90
            extensions/density_qn/tests/density_qn_mock.f90
            extensions/density_qn/tests/density_qn_c_interface_mock.f90
            extensions/density_qn/tests/density_qn_c_integration_tests.c
        )
    endif()
    if(ENABLE_S_GEK)
        list(APPEND OPENTRUSTREGION_TESTS
            extensions/s_gek/tests/s_gek_test_reference.f90
            extensions/s_gek/tests/s_gek_unit_tests.f90
            extensions/s_gek/tests/s_gek_c_interface_unit_tests.f90
            extensions/s_gek/tests/s_gek_mock.f90
            extensions/s_gek/tests/s_gek_c_interface_mock.f90
            extensions/s_gek/tests/s_gek_c_integration_tests.c
        )
    endif()

    # create shared test library which is then called from Python for testing
    add_library(otrtestsuite SHARED ${OPENTRUSTREGION_TESTS}
                $<TARGET_OBJECTS:opentrustregion>)

    # enable preprocessing
    set_source_files_properties(${OPENTRUSTREGION_TESTS}
                                PROPERTIES Fortran_PREPROCESS ON
    )

    # link against BLAS and LAPACK and against OpenTrustRegion
    target_link_libraries(otrtestsuite
                          PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} opentrustregion
    )

    # install testsuite
    if (${otr}_INSTALL_TESTING_LIBRARY)
        install(TARGETS otrtestsuite
                RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT ${test_component} # Windows (.dll)
                LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT ${test_component} # shared libraries (.so, .dylib)
                ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT ${test_component} # static libraries (.a)
        )
    endif()
endif()

# create targets file in the build tree for build-tree usage
export(TARGETS opentrustregion
       FILE "${CMAKE_CURRENT_BINARY_DIR}/${otr}Targets-${LIB_SUFFIX}.cmake"
       NAMESPACE "${otr}::"
)

# generate OpenTrustRegionConfig.cmake
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/${otr}Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${otr}Config.cmake
    INSTALL_DESTINATION ${${otr}_INSTALL_CMAKEDIR}
)

# generate version file, which accepts only the requested minor version, since a
# program has to be rebuilt for another one
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/${otr}ConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMinorVersion
)

# install library
install(TARGETS opentrustregion
        EXPORT library_interface
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT ${runtime_component} # Windows (.dll)
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT ${runtime_component} # shared libraries (.so, .dylib)
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT ${runtime_component} # static libraries (.a)
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT ${runtime_component}
)

# install headers
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        COMPONENT ${runtime_component}
)
if(ANY_EXTENSION_ENABLED)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/extensions/common/include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            COMPONENT ${runtime_component}
    )
endif()
if(ENABLE_MO)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/extensions/mo/include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            COMPONENT ${runtime_component}
    )
endif()
if(ENABLE_OAO)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/extensions/oao/include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            COMPONENT ${runtime_component}
    )
endif()
if(ENABLE_ROTATION_QN)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/extensions/rotation_qn/include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            COMPONENT ${runtime_component}
    )
endif()
if(ENABLE_DENSITY_QN)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/extensions/density_qn/include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            COMPONENT ${runtime_component}
    )
endif()
if(ENABLE_S_GEK)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/extensions/s_gek/include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            COMPONENT ${runtime_component}
    )
endif()

# install Fortran modules
install(DIRECTORY ${module_dir}/
        DESTINATION ${module_install_dir}
        COMPONENT ${runtime_component}
        FILES_MATCHING PATTERN "*.mod"
)

# install targets
install(EXPORT library_interface
        FILE ${otr}Targets-${LIB_SUFFIX}.cmake
        NAMESPACE "${otr}::"
        DESTINATION ${${otr}_INSTALL_CMAKEDIR}
        COMPONENT ${runtime_component}
)

# install config files
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${otr}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${otr}ConfigVersion.cmake
        DESTINATION ${${otr}_INSTALL_CMAKEDIR}
        COMPONENT ${runtime_component}
)

# ====  Python wheel  ==========================================================

# when the build is driven by scikit-build-core, install the compiled libraries
# directly into the Python package directory so that ctypes finds them next to
# python_interface.py
if(SKBUILD)
    # a wheel stores plain files, so unset the version properties to drop the versioned
    # symlink chain and emit a single library named exactly as the ctypes loader expects
    set_property(TARGET opentrustregion PROPERTY VERSION)
    set_property(TARGET opentrustregion PROPERTY SOVERSION)

    # keep the directories of the linked libraries outside the system's, such as a
    # BLAS/LAPACK installed elsewhere, in the installed libraries, as in the build
    # tree, so that a source install finds them without LD_LIBRARY_PATH; the wheels
    # link system libraries only, which CMake leaves out
    set_property(TARGET opentrustregion PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

    # install the solver library
    install(TARGETS opentrustregion
            RUNTIME DESTINATION ${SKBUILD_PLATLIB_DIR}/pyopentrustregion COMPONENT ${python_component} # Windows (.dll)
            LIBRARY DESTINATION ${SKBUILD_PLATLIB_DIR}/pyopentrustregion COMPONENT ${python_component} # shared libraries (.so, .dylib)
    )

    # install the testsuite library whenever it was built, so that an editable or
    # source install keeps "python -m pyopentrustregion.testsuite" working; release
    # wheels switch the testsuite off instead of excluding it here
    if(${otr}_BUILD_TESTING AND NOT ${otr}_HOST_PROVIDES_BLAS)
        set_property(TARGET otrtestsuite PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
        install(TARGETS otrtestsuite
                RUNTIME DESTINATION ${SKBUILD_PLATLIB_DIR}/pyopentrustregion COMPONENT ${python_component} # Windows (.dll)
                LIBRARY DESTINATION ${SKBUILD_PLATLIB_DIR}/pyopentrustregion COMPONENT ${python_component} # shared libraries (.so, .dylib)
        )
    endif()
endif()
