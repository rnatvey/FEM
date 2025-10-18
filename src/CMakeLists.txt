file(GLOB SRC_HEADERS include/*.h)
file(GLOB SRC_SOURCES src/*.cpp)


add_library(src)
target_sources(src PUBLIC 
    src/main.cpp
)



target_link_libraries(src 
linalg
solver
)
add_executable(FEM src/main.cpp)

target_include_directories(src PRIVATE
    ${CMAKE_SOURCE_DIR}/Libs/solver/include  
)
target_link_libraries(FEM solver linalg)

target_compile_features(solver PUBLIC cxx_std_23)
target_compile_definitions(solver PUBLIC SOLVER_EXPORT)


