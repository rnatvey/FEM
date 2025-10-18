cmake_minimum_required(VERSION 3.15)
project(FEM)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(Libs/linalg)
add_subdirectory(Libs/solver)
add_subdirectory(src)

if(BUIL_TESTS)
	add_subdirectory(tests)
endif()

if(MSVC)
    add_compile_options(/showIncludes-)
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY CLEAN_NO_CUSTOM 1)