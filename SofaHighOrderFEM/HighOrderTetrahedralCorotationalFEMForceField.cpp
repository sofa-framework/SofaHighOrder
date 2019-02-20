
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
        .add< HighOrderTetrahedralCorotationalFEMForceField<Vec3Types> >()
        ;

template class  SOFA_HIGHORDER_FEM_API HighOrderTetrahedralCorotationalFEMForceField<Vec3Types>;

} // namespace forcefield

} // namespace component

} // namespace sofa

