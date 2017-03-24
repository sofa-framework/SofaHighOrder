
#include <sofa/helper/system/config.h>
#include "initHighOrderTopology.h"

namespace sofa
{

namespace component
{
//Here are just several convenient functions to help users know what the plugin contains 

extern "C" {
    SOFA_HIGHORDER_TOPOLOGY_API void initExternalModule();
    SOFA_HIGHORDER_TOPOLOGY_API const char* getModuleName();
    SOFA_HIGHORDER_TOPOLOGY_API const char* getModuleVersion();
    SOFA_HIGHORDER_TOPOLOGY_API const char* getModuleLicense();
    SOFA_HIGHORDER_TOPOLOGY_API const char* getModuleDescription();
    SOFA_HIGHORDER_TOPOLOGY_API const char* getModuleComponentList();
}

void initExternalModule()
{
    static bool first = true;
    if (first)
    {
        first = false;
    }
}

const char* getModuleName()
{
    return "HighOrderTopology Plugin";
}

const char* getModuleVersion()
{
    return "0.5";
}

const char* getModuleLicense()
{
    return "LGPL";
}


const char* getModuleDescription()
{
    return "High Order Triangles, Tetrahedra Topology classes into SOFA";
}

const char* getModuleComponentList()
{
    return "HighOrderTetrahedronSetTopologyContainer, BezierTetrahedronSetGeometryAlgorithms, LagrangeTetrahedronSetGeometryAlgorithms,HighOrderTriangleSetTopologyContainer, BezierTriangleSetGeometryAlgorithms, LagrangeTriangleSetGeometryAlgorithms,GenerateBezierSphere, GenerateBezierCylinder";
}




} // namespace component

} // namespace sofa
