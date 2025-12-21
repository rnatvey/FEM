file(GLOB MESH_HEADERS include/*.h)
file(GLOB MESH_SOURCES src/*.cpp)


add_library(mesh)
target_sources(mesh PUBLIC  
    "assembly/include/assembly.h"
    "element/include/BaseElement.h"
    "element/include/node.h"
    "element/include/material.h"
    "element/include/planeisometric/ContactPlane.h"
    #"element/include/planeisometric/StaticPlane.h"
    "element/include/planeisometric/Planeisoparametric.h"
    "element/include/planeisometric/Planeisoparametric.cpp"
    "element/include/planeisometric/ContactPlane.cpp"
    "assembly/src/assembly.cpp"
    "element/src/BaseElement.cpp"
    "element/src/material.cpp"
    "model/include/FEMModel.h"
    "model/src/FEMModel.cpp"
    "model/include/ConcentratedForce.h" 
    "model/src/ConcentratedForce.cpp"
    "model/include/geometry.h"
    "meshGenerator/include/meshgenerator.h"
    "meshGenerator/src/meshgenerator.cpp"
    "model/src/geometry.cpp" 
    "model/src/loadFunction.cpp")

target_include_directories(mesh PUBLIC
 ${CMAKE_CURRENT_SOURCE_DIR}/element/include
 ${CMAKE_CURRENT_SOURCE_DIR}/meshing/include
 ${CMAKE_CURRENT_SOURCE_DIR}/assembly/include
 ${CMAKE_CURRENT_SOURCE_DIR}/model/include
 ${CMAKE_CURRENT_SOURCE_DIR}/meshGenerator/include
)

target_include_directories(solver PRIVATE
    ${CMAKE_SOURCE_DIR}/Libs/mesh/element/include 
    ${CMAKE_SOURCE_DIR}/Libs/mesh/assembly/include
    ${CMAKE_SOURCE_DIR}/Libs/mesh/meshing/include 
    ${CMAKE_SOURCE_DIR}/Libs/solver/include
    ${CMAKE_SOURCE_DIR}/3rdParty/eigen-5.0.0 
    ${CMAKE_SOURCE_DIR}/3rdParty/boost/libs/math/include
)

target_link_libraries(mesh
solver
Eigen3::Eigen
Boost::math
)

target_compile_features(mesh PUBLIC cxx_std_23)
target_compile_definitions(mesh PUBLIC SOLVER_EXPORT)