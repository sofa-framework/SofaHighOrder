cmake_minimum_required(VERSION 2.8)
project(SofaHighOrderFEM)

set(SOFAHIGHORDERFEM_MAJOR_VERSION 0)
set(SOFAHIGHORDERFEM_MINOR_VERSION 2)
set(SOFAHIGHORDERFEM_VERSION ${SOFAHIGHORDERFEM_MAJOR_VERSION}.${SOFAHIGHORDERFEM_MINOR_VERSION})

set(HEADER_FILES
	HighOrderMeshMatrixMass.h
	HighOrderMeshMatrixMass.inl
    initHighOrderFEM.h
	AdaptiveBezierMeshMatrixMass.h
	AdaptiveBezierMeshMatrixMass.inl
	AdaptiveBezierTetrahedronController.h
	AdaptiveBezierTetrahedronController.inl
	AdaptiveBezierTetra2BezierTriangleTopologicalMapping.h
	AdaptiveBezierTetra2BezierTriangleMechanicalMapping.h
	AdaptiveBezierTetra2BezierTriangleMechanicalMapping.inl
	AdaptiveBezierTetrahedronSetTopologyModifier.h
	AdaptiveBezierTetrahedronSetTopologyContainer.h
	AdaptiveBezierTetrahedronSetTopologyAlgorithms.h
	AdaptiveBezierTetrahedronSetTopologyAlgorithms.inl
	AdaptiveBezierTetrahedralCorotationalFEMForceField.h 
    AdaptiveBezierTetrahedralCorotationalFEMForceField.inl
	HighOrderTetrahedralCorotationalFEMForceField.h 
    HighOrderTetrahedralCorotationalFEMForceField.inl
	HighOrderTriangularCorotationalFEMForceField.h 
    HighOrderTriangularCorotationalFEMForceField.inl
	HighOrderTetrahedralDiffusionForceField.h 
    HighOrderTetrahedralDiffusionForceField.inl
	HighOrderTriangularDiffusionForceField.h 
    HighOrderTriangularDiffusionForceField.inl
	)

set(SOURCE_FILES
	HighOrderMeshMatrixMass.cpp
	initHighOrderFEM.cpp
	AdaptiveBezierMeshMatrixMass.cpp
	AdaptiveBezierTetrahedronController.cpp
	AdaptiveBezierTetra2BezierTriangleTopologicalMapping.cpp
	AdaptiveBezierTetra2BezierTriangleMechanicalMapping.cpp
	AdaptiveBezierTetrahedronSetTopologyContainer.cpp
	AdaptiveBezierTetrahedronSetTopologyModifier.cpp
	AdaptiveBezierTetrahedronSetTopologyAlgorithms.cpp
    AdaptiveBezierTetrahedralCorotationalFEMForceField.cpp 
    HighOrderTetrahedralCorotationalFEMForceField.cpp 
    HighOrderTriangularCorotationalFEMForceField.cpp 
	HighOrderTetrahedralDiffusionForceField.cpp
	HighOrderTriangularDiffusionForceField.cpp
	)

find_package(SofaHighOrderTopology)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} )


set(README_FILE "SofaHighOrderFEM.txt" )
set(SOFAHIGHORDERFEM_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})

set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-DSOFA_BUILD_SOFAHIGHORDERFEM")
target_link_libraries(${PROJECT_NAME} SofaHighOrderTopology SofaComponentBase SofaBaseTopology SofaUserInteraction SofaMiscForceField)

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<INSTALL_INTERFACE:include/${PROJECT_NAME}>")


sofa_create_package(SofaHighOrderFEM ${SOFAHIGHORDERFEM_VERSION} SofaHighOrderFEM SofaHighOrderFEM)
#
# Organize projects into sofa-asclepios Visual Studio folders
#
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
SET_PROPERTY(TARGET ${PROJECT_NAME}	PROPERTY FOLDER "SofaHighOrderFEM")