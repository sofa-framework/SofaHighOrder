
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
        .add< HighOrderTriangularDiffusionForceField<Vec2Types> >()
        .add< HighOrderTriangularDiffusionForceField<Vec1Types> >()
        ;


template class  SOFA_HIGHORDER_FEM_API HighOrderTriangularDiffusionForceField<Vec2Types>;
template class  SOFA_HIGHORDER_FEM_API HighOrderTriangularDiffusionForceField<Vec1Types>;


} // namespace forcefield

} // namespace component

} // namespace sofa

