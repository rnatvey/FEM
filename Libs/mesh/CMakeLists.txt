file(GLOB MESH_HEADERS include/*.h)
file(GLOB MESH_SOURCES src/*.cpp)


add_library(mesh)
target_sources(mesh PUBLIC  
    "element/include/BaseElement.h"
)

target_include_directories(mesh PUBLIC
 ${CMAKE_CURRENT_SOURCE_DIR}/element/include
 ${CMAKE_CURRENT_SOURCE_DIR}/meshing/include
)

target_include_directories(solver PRIVATE
    ${CMAKE_SOURCE_DIR}/Libs/mesh/element/include  
    ${CMAKE_SOURCE_DIR}/Libs/mesh/meshing/include 
)

target_link_libraries(mesh PUBLIC linalg)
target_link_libraries(mesh PUBLIC solver)
target_compile_features(mesh PUBLIC cxx_std_23)
target_compile_definitions(mesh PUBLIC SOLVER_EXPORT)