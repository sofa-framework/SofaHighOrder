#define SOFA_COMPONENT_INTERACTIONFORCEFIELD_HIGHORDERTRIANGULARCOROTATIONALFEMFORCEFIELD_CPP
#include "HighOrderTriangularCorotationalFEMForceField.inl"
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>
#include "initHighOrderFEM.h"

namespace sofa
{

namespace component
{

namespace forcefield
{




SOFA_DECL_CLASS(HighOrderTriangularCorotationalFEMForceField)

using namespace sofa::defaulttype;


// Register in the Factory
int HighOrderTriangularCorotationalFEMForceFieldClass = core::RegisterObject("Corotational High Order Triangular Mesh")
        .add< HighOrderTriangularCorotationalFEMForceField<Vec2Types> >()
        ;

template class  SOFA_HIGHORDER_FEM_API HighOrderTriangularCorotationalFEMForceField<Vec2Types>;


} // namespace forcefield

} // namespace component

} // namespace sofa

