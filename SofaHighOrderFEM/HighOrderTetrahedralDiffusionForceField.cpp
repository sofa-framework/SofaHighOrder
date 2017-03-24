
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
#ifndef SOFA_FLOAT
        .add< HighOrderTetrahedralDiffusionForceField<Vec3dTypes> >()
        .add< HighOrderTetrahedralDiffusionForceField<Vec2dTypes> >()
        .add< HighOrderTetrahedralDiffusionForceField<Vec1dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< HighOrderTetrahedralDiffusionForceField<Vec3fTypes> >()
        .add< HighOrderTetrahedralDiffusionForceField<Vec2fTypes> >()
        .add< HighOrderTetrahedralDiffusionForceField<Vec1fTypes> >()
#endif
        ;

#ifndef SOFA_FLOAT
template class  SOFA_HIGHORDER_FEM_API HighOrderTetrahedralDiffusionForceField<Vec3dTypes>;
template class  SOFA_HIGHORDER_FEM_API HighOrderTetrahedralDiffusionForceField<Vec2dTypes>;
template class  SOFA_HIGHORDER_FEM_API HighOrderTetrahedralDiffusionForceField<Vec1dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_HIGHORDER_FEM_API HighOrderTetrahedralDiffusionForceField<Vec3fTypes>;
template class SOFA_HIGHORDER_FEM_API HighOrderTetrahedralDiffusionForceField<Vec2fTypes>;
template class SOFA_HIGHORDER_FEM_API HighOrderTetrahedralDiffusionForceField<Vec1fTypes>;
#endif

} // namespace forcefield

} // namespace component

} // namespace sofa

