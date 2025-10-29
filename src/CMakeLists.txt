file(GLOB SRC_HEADERS include/*.h)
file(GLOB SRC_SOURCES src/*.cpp)

set(CMAKE_SHOW_INCLUDES OFF)
set(CMAKE_INCLUDE_CURRENT_DIR OFF)

add_library(src)
target_sources(src PUBLIC 
    src/main.cpp
)
find_package(Eigen3 5.0 REQUIRED NO_MODULE)
target_link_libraries(src 
linalg
solver
Eigen3::Eigen
)
add_executable(FEM src/main.cpp)

target_include_directories(src PRIVATE
    ${CMAKE_SOURCE_DIR}/Libs/solver/include
    ${CMAKE_SOURCE_DIR}/Libs/linalg/include 
    ${CMAKE_SOURCE_DIR}/3rdParty/eigen-5.0.0 
)
target_link_libraries(FEM solver linalg Eigen3::Eigen)

target_compile_features(solver PUBLIC cxx_std_23)
target_compile_definitions(solver PUBLIC SOLVER_EXPORT)


