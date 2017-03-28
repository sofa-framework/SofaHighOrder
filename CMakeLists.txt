cmake_minimum_required(VERSION 2.8.12)

find_package(SofaFramework REQUIRED)
sofa_add_plugin(SofaHighOrderTopology/ SofaHighOrderTopology)
sofa_add_plugin(SofaHighOrderFEM/ SofaHighOrderFEM)


sofa_install_git_version( "plugin.HighOrder" ${CMAKE_CURRENT_SOURCE_DIR} )
