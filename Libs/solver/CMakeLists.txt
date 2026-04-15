add_library(solver)

target_sources(solver PRIVATE
    src/solver.cpp
)

target_link_libraries(solver PUBLIC
    linalg
    Eigen3::Eigen
)

target_include_directories(solver PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if(USE_INTEL_MKL)
    target_link_libraries(solver PUBLIC IntelMKL)
    target_compile_definitions(solver PUBLIC SOLVER_USE_MKL)
endif()

target_compile_options(solver PRIVATE
    $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU>>:-O3 -march=native -mtune=native>
    $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:Clang>>:-O3 -march=native -mtune=native>
    $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:MSVC>>:/O2 /arch:AVX2>
)

option(ENABLE_EIGEN_DEBUG "Enable Eigen runtime checks" OFF)
if(NOT ENABLE_EIGEN_DEBUG)
    target_compile_definitions(solver PRIVATE EIGEN_NO_DEBUG)
endif()

target_compile_features(solver PUBLIC cxx_std_23)