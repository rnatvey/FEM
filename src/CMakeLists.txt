file(GLOB SRC_HEADERS include/*.h)
file(GLOB SRC_SOURCES src/*.cpp)

set(CMAKE_SHOW_INCLUDES OFF)
set(CMAKE_INCLUDE_CURRENT_DIR OFF)

add_executable(FEM 
    src/main.cpp
)

target_include_directories(FEM PRIVATE
    ${CMAKE_SOURCE_DIR}/Libs/solver/include
    ${CMAKE_SOURCE_DIR}/Libs/linalg/include 
    ${CMAKE_SOURCE_DIR}/Libs/mesh/element/include 
    ${CMAKE_SOURCE_DIR}/Libs/mesh/assembly/include 
    ${CMAKE_SOURCE_DIR}/Libs/mesh/meshGenerator/include 
    ${CMAKE_SOURCE_DIR}/3rdParty/eigen-5.0.0 
)

target_link_libraries(FEM 
    PRIVATE
        solver
        linalg
        mesh
        Eigen3::Eigen
)

if(USE_INTEL_MKL)
    target_link_libraries(FEM PRIVATE IntelMKL)
    target_compile_definitions(FEM PRIVATE
        EIGEN_USE_MKL_ALL
        EIGEN_USE_MKL
        SOLVER_USE_MKL
    )
    
    if(WIN32)
        # Получаем путь к MKL
        if(NOT DEFINED MKLROOT AND DEFINED ENV{MKLROOT})
            set(MKLROOT $ENV{MKLROOT})
        endif()
        
        if(DEFINED MKLROOT)
            # ИСПРАВЛЕННЫЕ ПУТИ с учетом вашей структуры
            set(MKL_DLL_SEARCH_PATHS
                "${MKLROOT}/redist/intel64"
                "${MKLROOT}/../redist/intel64"
                "C:/Program Files (x86)/Intel/oneAPI/redist/intel64"
                "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/redist/intel64"
            )
            
            # Ищем MKL DLL
            set(MKL_DLL_NAMES
                mkl_core.dll
                mkl_intel_thread.dll
                mkl_intel_lp64.dll
            )
            
            set(MKL_DLLS_FOUND)
            
            foreach(DLL_NAME ${MKL_DLL_NAMES})
                find_file(MKL_DLL_PATH_${DLL_NAME}
                    NAMES ${DLL_NAME}
                    PATHS ${MKL_DLL_SEARCH_PATHS}
                    NO_DEFAULT_PATH
                )
                
                if(MKL_DLL_PATH_${DLL_NAME})
                    list(APPEND MKL_DLLS_FOUND "${MKL_DLL_PATH_${DLL_NAME}}")
                    message(STATUS "Found ${DLL_NAME}: ${MKL_DLL_PATH_${DLL_NAME}}")
                else()
                    message(WARNING "Could not find ${DLL_NAME}")
                endif()
            endforeach()
            
            # ИСПРАВЛЕННЫЙ ПУТЬ для OpenMP DLL - используем вашу конкретную версию
            set(OMP_DLL_SEARCH_PATHS
                "C:/Program Files (x86)/Intel/oneAPI/compiler/2025.3/bin"
                "C:/Program Files (x86)/Intel/oneAPI/compiler/latest/bin"
                "${MKLROOT}/../compiler/2025.3/bin"
                "${MKLROOT}/../../compiler/2025.3/bin"
                "C:/Program Files (x86)/Intel/oneAPI/compiler/2025.3/lib/intel64"
                "C:/Program Files (x86)/Intel/oneAPI/compiler/latest/lib/intel64"
            )
            
            find_file(OMP_DLL_PATH
                NAMES libiomp5md.dll
                PATHS ${OMP_DLL_SEARCH_PATHS}
                NO_DEFAULT_PATH
            )
            
            if(OMP_DLL_PATH)
                list(APPEND MKL_DLLS_FOUND "${OMP_DLL_PATH}")
                message(STATUS "Found OpenMP DLL: ${OMP_DLL_PATH}")
            else()
                message(WARNING "Could not find libiomp5md.dll in ${OMP_DLL_SEARCH_PATHS}")
            endif()
            
            # Копируем все найденные DLL
            foreach(DLL_FILE ${MKL_DLLS_FOUND})
                if(EXISTS "${DLL_FILE}")
                    message(STATUS "Copying DLL: ${DLL_FILE}")
                    add_custom_command(TARGET FEM POST_BUILD
                        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                            "${DLL_FILE}"
                            $<TARGET_FILE_DIR:FEM>
                        COMMENT "Copying ${DLL_FILE} to output directory"
                    )
                endif()
            endforeach()
            
            # Если не нашли DLL, показываем где искали
            if(NOT MKL_DLLS_FOUND)
                message(WARNING "No MKL DLLs found. Searched in: ${MKL_DLL_SEARCH_PATHS}")
            endif()
        endif()
    endif()
endif()

target_compile_features(FEM PUBLIC cxx_std_23)