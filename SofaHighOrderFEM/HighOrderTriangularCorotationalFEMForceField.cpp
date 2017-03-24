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
#ifndef SOFA_FLOAT
        .add< HighOrderTriangularCorotationalFEMForceField<Vec2dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< HighOrderTriangularCorotationalFEMForceField<Vec2fTypes> >()
#endif
        ;

#ifndef SOFA_FLOAT
template class  SOFA_HIGHORDER_FEM_API HighOrderTriangularCorotationalFEMForceField<Vec2dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_HIGHORDER_FEM_API HighOrderTriangularCorotationalFEMForceField<Vec2fTypes>;
#endif

} // namespace forcefield

} // namespace component

} // namespace sofa

