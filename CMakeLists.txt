cmake_minimum_required(VERSION 3.15)
set(CMAKE_SHOW_INCLUDES OFF)
set(CMAKE_INCLUDE_CURRENT_DIR OFF)

project(FEM)

if(MSVC)
    set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY COMPILE_OPTIONS "/showIncludes")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(USE_INTEL_MKL "Use Intel MKL" ON)

# НАСТРОЙКА MKL - делаем ДО добавления поддиректорий
if(USE_INTEL_MKL)
    # Находим OpenMP
    find_package(OpenMP REQUIRED)
    message(STATUS "OpenMP found: ${OpenMP_CXX_FLAGS}")
    
    # Определяем путь к MKL
    if(NOT DEFINED ENV{MKLROOT} AND NOT DEFINED MKLROOT)
        set(MKLROOT "C:/Program Files (x86)/Intel/oneAPI/mkl/latest" CACHE PATH "Path to MKL")
    elseif(DEFINED ENV{MKLROOT})
        set(MKLROOT $ENV{MKLROOT} CACHE PATH "Path to MKL")
    endif()
    
    if(EXISTS "${MKLROOT}/include/mkl.h")
        message(STATUS "MKL found at: ${MKLROOT}")
        
        # Библиотеки находятся прямо в lib (без intel64)
        set(MKL_LIB_PATH "${MKLROOT}/lib")
        
        # Проверяем существование библиотек
        if(EXISTS "${MKL_LIB_PATH}/mkl_intel_lp64.lib" AND 
           EXISTS "${MKL_LIB_PATH}/mkl_intel_thread.lib" AND 
           EXISTS "${MKL_LIB_PATH}/mkl_core.lib")
            
            message(STATUS "MKL libraries found in: ${MKL_LIB_PATH}")
            
            # Создаем интерфейсную библиотеку для MKL
            add_library(IntelMKL INTERFACE)
            
            # Добавляем include директории
            target_include_directories(IntelMKL INTERFACE 
                ${MKLROOT}/include
            )
            
            # Добавляем библиотеки для линковки (полные пути)
            target_link_libraries(IntelMKL INTERFACE 
                "${MKL_LIB_PATH}/mkl_intel_lp64.lib"
                "${MKL_LIB_PATH}/mkl_intel_thread.lib"
                "${MKL_LIB_PATH}/mkl_core.lib"
            )
            
            # Ищем OpenMP библиотеку Intel
            set(INTEL_OMP_PATH "${MKLROOT}/../compiler/latest/lib")
            if(NOT EXISTS "${INTEL_OMP_PATH}/libiomp5md.lib")
                set(INTEL_OMP_PATH "${MKLROOT}/../../compiler/latest/lib")
            endif()
            
            if(EXISTS "${INTEL_OMP_PATH}/libiomp5md.lib")
                message(STATUS "Intel OMP library found at: ${INTEL_OMP_PATH}/libiomp5md.lib")
                target_link_libraries(IntelMKL INTERFACE 
                    "${INTEL_OMP_PATH}/libiomp5md.lib"
                )
            else()
                message(WARNING "Intel OMP library not found. Using OpenMP::OpenMP_CXX")
                target_link_libraries(IntelMKL INTERFACE OpenMP::OpenMP_CXX)
            endif()
            
            # Добавляем определения компилятора
            target_compile_definitions(IntelMKL INTERFACE
                EIGEN_USE_MKL_ALL
                EIGEN_USE_MKL
                NOMINMAX
                _CRT_SECURE_NO_WARNINGS
            )
            
            # Для MSVC добавляем флаги оптимизации
            target_compile_options(IntelMKL INTERFACE
                $<$<CXX_COMPILER_ID:MSVC>:/fp:fast /arch:AVX2 /openmp>
            )
            
            message(STATUS "IntelMKL interface library created successfully")
        else()
            message(FATAL_ERROR "MKL libraries not found in ${MKL_LIB_PATH}")
        endif()
    else()
        message(FATAL_ERROR "MKL not found at ${MKLROOT}")
    endif()
endif()

# ТЕПЕРЬ добавляем поддиректории
add_subdirectory(3rdParty/eigen-5.0.0)
add_subdirectory(Libs/linalg)
add_subdirectory(Libs/solver)
add_subdirectory(Libs/mesh)
add_subdirectory(src)

message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")

if(BUILD_TESTS)
    add_subdirectory(tests)
endif()

if(MSVC)
    add_compile_options(/showIncludes-)
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY CLEAN_NO_CUSTOM 1)