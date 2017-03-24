cmake_minimum_required(VERSION 2.8.12)

project(SofaHighOrderTopology_test)

set(SOURCE_FILES
    BezierTetrahedronTopology_test.cpp)

add_definitions("-DSOFAHIGHORDERTOPOLOGY_TEST_SCENES_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/scenes\"")
add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")
target_link_libraries(${PROJECT_NAME} SofaHighOrderTopology SofaGTestMain SofaTest)

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
