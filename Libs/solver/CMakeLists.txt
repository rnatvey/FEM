file(GLOB SOLVER_HEADERS include/*.h)
file(GLOB SOLVER_SOURCES src/*.cpp)

set(CMAKE_SHOW_INCLUDES OFF)
set(CMAKE_INCLUDE_CURRENT_DIR OFF)



add_library(solver)
target_sources(solver PUBLIC  
    "src/solver.cpp"
    "include/solver.h"
)

target_link_libraries(solver PUBLIC
linalg
Eigen3::Eigen
)

target_include_directories(solver PUBLIC
 ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_include_directories(solver PRIVATE
    ${CMAKE_SOURCE_DIR}/Libs/solver/include  
    ${CMAKE_SOURCE_DIR}/3rdParty/eigen-5.0.0 
)

if(USE_INTEL_MKL)
    target_link_libraries(solver PUBLIC IntelMKL)  # PUBLIC, чтобы зависимость передавалась дальше
    target_compile_definitions(solver PUBLIC SOLVER_USE_MKL)
endif()


set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")

# Оптимизации для Eigen
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(solver PRIVATE 
        $<$<CXX_COMPILER_ID:GNU>:-O3 -march=native -mtune=native>
        $<$<CXX_COMPILER_ID:Clang>:-O3 -march=native -mtune=native>
        $<$<CXX_COMPILER_ID:MSVC>:/O2 /arch:AVX2>
    )
endif()

# Для отладки (можно добавить опцию)
option(ENABLE_EIGEN_DEBUG "Enable Eigen debugging" OFF)
if(ENABLE_EIGEN_DEBUG)
    target_compile_definitions(solver PRIVATE EIGEN_NO_DEBUG)
endif()
target_compile_features(solver PUBLIC cxx_std_23)
target_compile_definitions(solver PUBLIC SOLVER_EXPORT)