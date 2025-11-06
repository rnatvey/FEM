file(GLOB SOLVER_HEADERS include/*.h)
file(GLOB SOLVER_SOURCES src/*.cpp)

set(CMAKE_SHOW_INCLUDES OFF)
set(CMAKE_INCLUDE_CURRENT_DIR OFF)



add_library(solver)
target_sources(solver PUBLIC  
    "src/solver.cpp"
    "include/solver.h"
)

target_link_libraries(solver 
linalg
Eigen3::Eigen
Boost::math
)

target_include_directories(solver PUBLIC
 ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_include_directories(solver PRIVATE
    ${CMAKE_SOURCE_DIR}/Libs/solver/include  
    ${CMAKE_SOURCE_DIR}/3rdParty/eigen-5.0.0 
    ${CMAKE_SOURCE_DIR}/3rdParty/boost/libs/math/include
)

target_compile_features(solver PUBLIC cxx_std_23)
target_compile_definitions(solver PUBLIC SOLVER_EXPORT)