file(GLOB SOLVER_HEADERS include/*.h)
file(GLOB SOLVER_SOURCES src/*.cpp)

set(CMAKE_SHOW_INCLUDES OFF)
set(CMAKE_INCLUDE_CURRENT_DIR OFF)

add_library(solver)
target_sources(solver PUBLIC  
    src/slau.cpp
)

target_include_directories(solver PUBLIC
 ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_include_directories(solver PRIVATE
    ${CMAKE_SOURCE_DIR}/Libs/solver/include  
)

target_link_libraries(solver PUBLIC linalg)
target_compile_features(solver PUBLIC cxx_std_23)
target_compile_definitions(solver PUBLIC SOLVER_EXPORT)