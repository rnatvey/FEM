cmake_minimum_required(VERSION 3.15)
project(FEM)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories("Libs")
add_library( solver_lib "Libs/solver/src/slau.cpp")
add_executable(FEM "main.cpp" )
target_link_libraries(FEM solver_lib)