
#define SOFA_COMPONENT_INTERACTIONFORCEFIELD_BEZIERTETRAHEDRALCOROTATIONALFEMFORCEFIELD_CPP
#include "BezierTetrahedralCorotationalFEMForceField.inl"
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/Vec3Types.h>
#include "initHighOrderFEM.h"

namespace sofa
{

namespace component
{

namespace forcefield
{




SOFA_DECL_CLASS(BezierTetrahedralCorotationalFEMForceField)

using namespace sofa::defaulttype;


// Register in the Factory
int BezierTetrahedralCorotationalFEMForceFieldClass = core::RegisterObject("Corotational Bezier Tetrahedral Mesh")
#ifndef SOFA_FLOAT
        .add< BezierTetrahedralCorotationalFEMForceField<Vec3dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< BezierTetrahedralCorotationalFEMForceField<Vec3fTypes> >()
#endif
        ;

#ifndef SOFA_FLOAT
template class  SOFA_HIGHORDER_FEM_API BezierTetrahedralCorotationalFEMForceField<Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_HIGHORDER_FEM_API BezierTetrahedralCorotationalFEMForceField<Vec3fTypes>;
#endif

} // namespace forcefield

} // namespace component

} // namespace sofa

