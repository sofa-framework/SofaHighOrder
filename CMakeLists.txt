cmake_minimum_required(VERSION 2.8.12)

find_package(SofaFramework REQUIRED)
sofa_add_plugin(SofaHighOrderTopology/ SofaHighOrderTopology)
sofa_add_plugin(SofaHighOrderFEM/ SofaHighOrderFEM)

if(SOFA_BUILD_TESTS)
    find_package(SofaTest QUIET)
    if(SofaTest_FOUND)
        add_subdirectory(SofaHighOrderTopology_test)
    endif()
endif()
#sofa_install_git_version( "sofa-dev" ${CMAKE_CURRENT_SOURCE_DIR} )
