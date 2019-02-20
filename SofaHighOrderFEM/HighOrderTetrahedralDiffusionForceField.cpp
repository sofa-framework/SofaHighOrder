
#define SOFA_COMPONENT_INTERACTIONFORCEFIELD_HIGHORDERTETRAHEDRALCOROTATIONALFEMFORCEFIELD_CPP
#include "HighOrderTetrahedralDiffusionForceField.inl"
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/Vec3Types.h>
#include "initHighOrderFEM.h"

namespace sofa
{

namespace component
{

namespace forcefield
{




SOFA_DECL_CLASS(HighOrderTetrahedralDiffusionForceField)

using namespace sofa::defaulttype;


// Register in the Factory
int HighOrderTetrahedralDiffusionForceFieldClass = core::RegisterObject("Diffusion on High Order Tetrahedral Mesh")
        .add< HighOrderTetrahedralDiffusionForceField<Vec3Types> >()
        .add< HighOrderTetrahedralDiffusionForceField<Vec2Types> >()
        .add< HighOrderTetrahedralDiffusionForceField<Vec1Types> >()
        ;


template class  SOFA_HIGHORDER_FEM_API HighOrderTetrahedralDiffusionForceField<Vec3Types>;
template class  SOFA_HIGHORDER_FEM_API HighOrderTetrahedralDiffusionForceField<Vec2Types>;
template class  SOFA_HIGHORDER_FEM_API HighOrderTetrahedralDiffusionForceField<Vec1Types>;


} // namespace forcefield

} // namespace component

} // namespace sofa

