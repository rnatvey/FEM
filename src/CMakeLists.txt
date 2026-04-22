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
    target_compile_definitions(${target} PRIVATE FEM_REPO_ROOT="${CMAKE_SOURCE_DIR}")

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

fem_configure_app(FEMBasicLinearReference src/basic_linear_reference.cpp)
fem_configure_app(FEMHyperelasticNoContactSanity src/hyperelastic_no_contact_sanity.cpp)
fem_configure_app(FEMHyperelasticLargeDeformationNoContact src/hyperelastic_large_deformation_no_contact.cpp)
fem_configure_app(FEMHyperelasticBlockOnRigidPlane src/hyperelastic_block_on_rigid_plane.cpp)
fem_configure_app(FEMHyperelasticContactGuard src/hyperelastic_contact_guard.cpp)
fem_configure_app(FEMSmallStrainStressSemanticsRegression src/small_strain_stress_semantics_regression.cpp)
fem_configure_app(FEM src/tire_contact_single_case.cpp)
fem_configure_app(FEMContactFarPlaneScenario src/contact_far_plane_scenario.cpp)
fem_configure_app(FEMBlockOnRigidPlane src/block_on_rigid_plane.cpp)
fem_configure_app(FEMContactBlockAugmentedLagrangianRegression src/contact_block_augmented_lagrangian_regression.cpp)
fem_configure_app(FEMContactBlockRegression src/contact_block_regression.cpp)
fem_configure_app(FEMContactNoContactAugmentedLagrangianRegression src/contact_no_contact_augmented_lagrangian_regression.cpp)
fem_configure_app(FEMContactNoContactRegression src/contact_no_contact_regression.cpp)
fem_configure_app(FEMRingContactStudy src/ring_contact_study.cpp)
fem_configure_app(FEMMainScaleContactStudy src/main_scale_contact_study.cpp)
fem_configure_app(FEMMainScaleAugmentedLagrangianStudy src/main_scale_augmented_lagrangian_study.cpp)
fem_configure_app(FEMMainScaleContactMethodComparisonStudy src/main_scale_contact_method_comparison_study.cpp)

enable_testing()
add_test(NAME FEMBasicLinearReference COMMAND FEMBasicLinearReference)
add_test(NAME FEMHyperelasticNoContactSanity COMMAND FEMHyperelasticNoContactSanity)
add_test(NAME FEMHyperelasticLargeDeformationNoContact COMMAND FEMHyperelasticLargeDeformationNoContact)
add_test(NAME FEMHyperelasticBlockOnRigidPlane COMMAND FEMHyperelasticBlockOnRigidPlane)
add_test(NAME FEMHyperelasticContactGuard COMMAND FEMHyperelasticContactGuard)
add_test(NAME FEMSmallStrainStressSemanticsRegression COMMAND FEMSmallStrainStressSemanticsRegression)
add_test(NAME FEMContactFarPlaneScenario COMMAND FEMContactFarPlaneScenario)
add_test(NAME FEMContactBlockAugmentedLagrangianRegression COMMAND FEMContactBlockAugmentedLagrangianRegression)
add_test(NAME FEMContactBlockRegression COMMAND FEMContactBlockRegression)
add_test(NAME FEMContactNoContactAugmentedLagrangianRegression COMMAND FEMContactNoContactAugmentedLagrangianRegression)
add_test(NAME FEMContactNoContactRegression COMMAND FEMContactNoContactRegression)
