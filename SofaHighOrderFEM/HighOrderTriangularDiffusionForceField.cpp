
#define SOFA_COMPONENT_INTERACTIONFORCEFIELD_HIGHORDERTRIANGULARCOROTATIONALFEMFORCEFIELD_CPP
#include "HighOrderTriangularDiffusionForceField.inl"
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>
#include "initHighOrderFEM.h"

namespace sofa
{

namespace component
{

namespace forcefield
{




SOFA_DECL_CLASS(HighOrderTriangularDiffusionForceField)

using namespace sofa::defaulttype;


// Register in the Factory
int HighOrderTriangularDiffusionForceFieldClass = core::RegisterObject("Diffusion on High Order Triangular Mesh")
#ifndef SOFA_FLOAT
        .add< HighOrderTriangularDiffusionForceField<Vec2dTypes> >()
        .add< HighOrderTriangularDiffusionForceField<Vec1dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< HighOrderTriangularDiffusionForceField<Vec2fTypes> >()
        .add< HighOrderTriangularDiffusionForceField<Vec1fTypes> >()
#endif
        ;

#ifndef SOFA_FLOAT
template class  SOFA_HIGHORDER_FEM_API HighOrderTriangularDiffusionForceField<Vec2dTypes>;
template class  SOFA_HIGHORDER_FEM_API HighOrderTriangularDiffusionForceField<Vec1dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_HIGHORDER_FEM_API HighOrderTriangularDiffusionForceField<Vec2fTypes>;
template class SOFA_HIGHORDER_FEM_API HighOrderTriangularDiffusionForceField<Vec1fTypes>;
#endif

} // namespace forcefield

} // namespace component

} // namespace sofa

