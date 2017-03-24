#define SOFA_COMPONENT_INTERACTIONFORCEFIELD_ADAPTIVEBEZIERTETRAHEDRALCOROTATIONALFEMFORCEFIELD_CPP
#include "AdaptiveBezierTetrahedralCorotationalFEMForceField.inl"
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/Vec3Types.h>
#include "initHighOrderFEM.h"

namespace sofa
{

namespace component
{

namespace forcefield
{

using namespace sofa::defaulttype;


SOFA_DECL_CLASS(AdaptiveBezierTetrahedralCorotationalFEMForceField)

using namespace sofa::defaulttype;


// Register in the Factory
int AdaptiveBezierTetrahedralCorotationalFEMForceFieldClass = core::RegisterObject("Adaptive Corotational Bezier Tetrahedral Mesh")
#ifndef SOFA_FLOAT
        .add< AdaptiveBezierTetrahedralCorotationalFEMForceField<Vec3dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< AdaptiveBezierTetrahedralCorotationalFEMForceField<Vec3fTypes> >()
#endif
        ;

#ifndef SOFA_FLOAT
template class  SOFA_HIGHORDER_FEM_API AdaptiveBezierTetrahedralCorotationalFEMForceField<Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_HIGHORDER_FEM_API AdaptiveBezierTetrahedralCorotationalFEMForceField<Vec3fTypes>;
#endif

} // namespace forcefield

} // namespace component

} // namespace sofa

