add_library(mesh)

target_sources(mesh PUBLIC
    assembly/include/assembly.h
    assembly/src/assembly.cpp
    element/include/BaseElement.h
    element/include/material.h
    element/include/node.h
    element/include/planeisometric/ContactPlane.h
    element/include/planeisometric/ContactPlane.cpp
    element/include/planeisometric/Planeisoparametric.h
    element/include/planeisometric/Planeisoparametric.cpp
    element/src/BaseElement.cpp
    element/src/material.cpp
    meshGenerator/include/meshgenerator.h
    meshGenerator/src/meshgenerator.cpp
    model/include/ConcentratedForce.h
    model/include/constants.h
    model/include/ContactTypes.h
    model/include/FEMModel.h
    model/include/geometry.h
    model/include/RigidPlaneContactSolver.h
    model/src/ConcentratedForce.cpp
    model/src/FEMModel.cpp
    model/src/geometry.cpp
    model/src/loadFunction.cpp
    model/src/RigidPlaneContactSolver.cpp
)

target_include_directories(mesh PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/assembly/include
    ${CMAKE_CURRENT_SOURCE_DIR}/element/include
    ${CMAKE_CURRENT_SOURCE_DIR}/meshGenerator/include
    ${CMAKE_CURRENT_SOURCE_DIR}/model/include
)

if(USE_INTEL_MKL)
    target_link_libraries(mesh PUBLIC IntelMKL)
    target_compile_definitions(mesh PUBLIC SOLVER_USE_MKL)
endif()

target_link_libraries(mesh PUBLIC
    solver
    Eigen3::Eigen
)

target_compile_features(mesh PUBLIC cxx_std_23)
target_compile_definitions(mesh PUBLIC SOLVER_EXPORT)
