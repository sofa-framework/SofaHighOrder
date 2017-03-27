#include <sofa/helper/system/config.h>
#include "initHighOrderFEM.h"


namespace sofa
{

namespace component
{


	extern "C" {
		SOFA_HIGHORDER_FEM_API void initExternalModule();
		SOFA_HIGHORDER_FEM_API const char* getModuleName();
		SOFA_HIGHORDER_FEM_API const char* getModuleVersion();
		SOFA_HIGHORDER_FEM_API const char* getModuleLicense();
		SOFA_HIGHORDER_FEM_API const char* getModuleDescription();
		SOFA_HIGHORDER_FEM_API const char* getModuleComponentList();
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
    return "HighOrderFEM Plugin";
}

const char* getModuleVersion()
{
    return "0.5";
}

const char* getModuleLicense()
{
    return "private";
}


const char* getModuleDescription()
{
    return "Elasticity on High Order Triangles, Tetrahedra";
}

const char* getModuleComponentList()
{
    return "HighOrderMeshMatrixMass, HighOrderTetrahedralCorotationalFEMForceField";
}



SOFA_LINK_CLASS(HighOrderTetrahedralCorotationalFEMForceField)
SOFA_LINK_CLASS(HighOrderMeshMatrixMass)

} // namespace component

} // namespace sofa
