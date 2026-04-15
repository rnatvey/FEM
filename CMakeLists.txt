cmake_minimum_required(VERSION 3.15)

project(FEM LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(USE_INTEL_MKL "Use Intel MKL-backed Eigen/Pardiso solver support" ON)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if(USE_INTEL_MKL)
    find_package(OpenMP REQUIRED)

    if(NOT DEFINED MKLROOT)
        if(DEFINED ENV{MKLROOT})
            set(MKLROOT "$ENV{MKLROOT}" CACHE PATH "Path to Intel oneAPI MKL")
        else()
            set(MKLROOT "C:/Program Files (x86)/Intel/oneAPI/mkl/latest" CACHE PATH "Path to Intel oneAPI MKL")
        endif()
    endif()

    set(MKL_INCLUDE_DIR "${MKLROOT}/include")
    set(MKL_LIBRARY_DIR "${MKLROOT}/lib")
    set(MKL_REQUIRED_LIBS
        "${MKL_LIBRARY_DIR}/mkl_intel_lp64.lib"
        "${MKL_LIBRARY_DIR}/mkl_intel_thread.lib"
        "${MKL_LIBRARY_DIR}/mkl_core.lib"
    )

    if(NOT EXISTS "${MKL_INCLUDE_DIR}/mkl.h")
        message(FATAL_ERROR "MKL headers were not found at ${MKL_INCLUDE_DIR}")
    endif()

    foreach(_mkl_lib IN LISTS MKL_REQUIRED_LIBS)
        if(NOT EXISTS "${_mkl_lib}")
            message(FATAL_ERROR "Required MKL library was not found: ${_mkl_lib}")
        endif()
    endforeach()

    set(INTEL_OMP_LIBRARY_CANDIDATES
        "${MKLROOT}/../compiler/latest/lib/libiomp5md.lib"
        "${MKLROOT}/../../compiler/latest/lib/libiomp5md.lib"
        "C:/Program Files (x86)/Intel/oneAPI/compiler/latest/lib/libiomp5md.lib"
    )

    set(INTEL_OMP_LIBRARY "")
    foreach(_omp_lib IN LISTS INTEL_OMP_LIBRARY_CANDIDATES)
        if(EXISTS "${_omp_lib}")
            set(INTEL_OMP_LIBRARY "${_omp_lib}")
            break()
        endif()
    endforeach()

    add_library(IntelMKL INTERFACE)
    target_include_directories(IntelMKL INTERFACE ${MKL_INCLUDE_DIR})
    target_link_libraries(IntelMKL INTERFACE ${MKL_REQUIRED_LIBS})

    if(INTEL_OMP_LIBRARY)
        target_link_libraries(IntelMKL INTERFACE "${INTEL_OMP_LIBRARY}")
    else()
        target_link_libraries(IntelMKL INTERFACE OpenMP::OpenMP_CXX)
    endif()

    target_compile_definitions(IntelMKL INTERFACE
        EIGEN_USE_MKL_ALL
        EIGEN_USE_MKL
        NOMINMAX
        _CRT_SECURE_NO_WARNINGS
    )

    target_compile_options(IntelMKL INTERFACE
        $<$<CXX_COMPILER_ID:MSVC>:/fp:fast /arch:AVX2 /openmp>
    )

    message(STATUS "MKL support enabled from ${MKLROOT}")
endif()

set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_DEMOS OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_CMAKE_PACKAGE OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_BLAS OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_LAPACK OFF CACHE BOOL "" FORCE)

add_subdirectory(3rdParty/eigen-5.0.0)
add_subdirectory(Libs/linalg)
add_subdirectory(Libs/solver)
add_subdirectory(Libs/mesh)
add_subdirectory(src)

message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")