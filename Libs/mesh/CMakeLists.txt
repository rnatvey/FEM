add_library(mesh)

target_sources(mesh PRIVATE
    assembly/src/assembly.cpp
    element/include/planeisometric/ContactPlane.cpp
    element/src/planeisometric/FiniteStrainQ4Element.cpp
    element/include/planeisometric/Planeisoparametric.cpp
    element/src/BaseElement.cpp
    element/src/NeoHookeanMaterial.cpp
    element/src/material.cpp
    meshGenerator/src/meshgenerator.cpp
    model/src/ConcentratedForce.cpp
    model/src/FEMModel.cpp
    model/src/geometry.cpp
    model/src/loadFunction.cpp
    model/src/ResultFileExporter.cpp
    model/src/RigidPlaneAugmentedLagrangianContactSolver.cpp
    model/src/RigidPlaneContactSolver.cpp
)

target_include_directories(mesh PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/assembly/include
    ${CMAKE_CURRENT_SOURCE_DIR}/element/include
    ${CMAKE_CURRENT_SOURCE_DIR}/meshGenerator/include
    ${CMAKE_CURRENT_SOURCE_DIR}/model/include
)

target_link_libraries(mesh PUBLIC
    solver
    Eigen3::Eigen
)

target_compile_features(mesh PUBLIC cxx_std_23)
