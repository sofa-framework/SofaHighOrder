
#define SOFA_COMPONENT_INTERACTIONFORCEFIELD_HIGHORDERTETRAHEDRALCOROTATIONALFEMFORCEFIELD_CPP
#include "HighOrderTetrahedralCorotationalFEMForceField.inl"
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/Vec3Types.h>
#include "initHighOrderFEM.h"

namespace sofa
{

namespace component
{

namespace forcefield
{




SOFA_DECL_CLASS(HighOrderTetrahedralCorotationalFEMForceField)

using namespace sofa::defaulttype;


// Register in the Factory
int HighOrderTetrahedralCorotationalFEMForceFieldClass = core::RegisterObject("Corotational High Order Tetrahedral Mesh")
#ifndef SOFA_FLOAT
        .add< HighOrderTetrahedralCorotationalFEMForceField<Vec3dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< HighOrderTetrahedralCorotationalFEMForceField<Vec3fTypes> >()
#endif
        ;

#ifndef SOFA_FLOAT
template class  SOFA_HIGHORDER_FEM_API HighOrderTetrahedralCorotationalFEMForceField<Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_HIGHORDER_FEM_API HighOrderTetrahedralCorotationalFEMForceField<Vec3fTypes>;
#endif

} // namespace forcefield

} // namespace component

} // namespace sofa

