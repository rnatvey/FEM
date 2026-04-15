add_library(linalg STATIC
    src/vector.cpp
)

target_include_directories(linalg PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(linalg PUBLIC cxx_std_23)