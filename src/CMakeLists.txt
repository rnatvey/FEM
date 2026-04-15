function(fem_copy_runtime_dll_if_exists target dll_path)
    if(EXISTS "${dll_path}")
        add_custom_command(TARGET ${target} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${dll_path}"
                $<TARGET_FILE_DIR:${target}>
            COMMENT "Copying ${dll_path} to output directory"
        )
    endif()
endfunction()

function(fem_configure_app target source_file)
    add_executable(${target} ${source_file})
    target_link_libraries(${target} PRIVATE mesh)
    target_compile_features(${target} PUBLIC cxx_std_23)

    if(USE_INTEL_MKL AND WIN32)
        set(_omp_dll_candidates
            "C:/Program Files (x86)/Intel/oneAPI/compiler/2025.3/bin/libiomp5md.dll"
            "C:/Program Files (x86)/Intel/oneAPI/compiler/latest/bin/libiomp5md.dll"
            "${MKLROOT}/../compiler/2025.3/bin/libiomp5md.dll"
            "${MKLROOT}/../../compiler/2025.3/bin/libiomp5md.dll"
            "${MKLROOT}/../compiler/latest/bin/libiomp5md.dll"
            "${MKLROOT}/../../compiler/latest/bin/libiomp5md.dll"
        )

        foreach(_dll_path IN LISTS _omp_dll_candidates)
            if(EXISTS "${_dll_path}")
                fem_copy_runtime_dll_if_exists(${target} "${_dll_path}")
                break()
            endif()
        endforeach()
    endif()
endfunction()

fem_configure_app(FEM src/main.cpp)
fem_configure_app(FEMBlockOnRigidPlane src/block_on_rigid_plane.cpp)
fem_configure_app(FEMContactBlockRegression src/contact_block_regression.cpp)
fem_configure_app(FEMContactNoContactRegression src/contact_no_contact_regression.cpp)
fem_configure_app(FEMRingContactStudy src/ring_contact_study.cpp)

enable_testing()
add_test(NAME FEMContactBlockRegression COMMAND FEMContactBlockRegression)
add_test(NAME FEMContactNoContactRegression COMMAND FEMContactNoContactRegression)
