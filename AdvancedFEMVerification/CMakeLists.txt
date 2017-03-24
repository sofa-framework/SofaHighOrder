cmake_minimum_required(VERSION 2.8)

project("AdvancedFEMVerification")

include(${SOFA_CMAKE_DIR}/preProject.cmake)

set(SOURCE_FILES
	AdvancedFEMVerification.cpp
	)



 
add_executable(${PROJECT_NAME} ${HEADER_FILES} ${SOURCE_FILES} ${RC_FILES})

AddCompilerDefinitions("SOFA_ADVANCEDFEM_VERIFICATION_SCENES_PATH=${CMAKE_CURRENT_SOURCE_DIR}")
AddLinkerDependencies(SofaAdvancedFEM 
                      SofaComponentMain
					  SofaSimulationTree
					  SofaSimulationGraph
					  )


include(${SOFA_CMAKE_DIR}/postProject.cmake)
