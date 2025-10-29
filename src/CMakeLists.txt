file(GLOB SRC_HEADERS include/*.h)
file(GLOB SRC_SOURCES src/*.cpp)

set(CMAKE_SHOW_INCLUDES OFF)
set(CMAKE_INCLUDE_CURRENT_DIR OFF)

add_library(src)
target_sources(src PUBLIC 
    src/main.cpp
)

target_link_libraries(src 
linalg
solver
Eigen3::Eigen
Boost::math
)
add_executable(FEM src/main.cpp)

target_include_directories(src PRIVATE
    ${CMAKE_SOURCE_DIR}/Libs/solver/include
    ${CMAKE_SOURCE_DIR}/Libs/linalg/include 
    ${CMAKE_SOURCE_DIR}/3rdParty/eigen-5.0.0 
    ${CMAKE_SOURCE_DIR}/3rdParty/boost/libs/math/include
)
target_link_libraries(FEM solver linalg Eigen3::Eigen Boost::math)
target_compile_features(solver PUBLIC cxx_std_23)
target_compile_definitions(solver PUBLIC SOLVER_EXPORT)


