#define SOFA_HIGHORDERTOPOLOGY_GENERATEBEZIERCYLINDER_CPP
#include "GenerateBezierCylinder.inl"
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/Vec3Types.h>

namespace sofa
{

namespace component
{

namespace engine
{
using namespace sofa::defaulttype;

SOFA_DECL_CLASS(GenerateBezierCylinder)

int GenerateBezierCylinderClass = core::RegisterObject("Generate a Cylindrical Tetrahedral Mesh")
#ifndef SOFA_FLOAT
        .add< GenerateBezierCylinder<Vec3dTypes> >()
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
        .add< GenerateBezierCylinder<Vec3fTypes> >()
#endif //SOFA_DOUBLE
        ;


#ifndef SOFA_FLOAT
template class SOFA_HIGHORDER_TOPOLOGY_API GenerateBezierCylinder<Vec3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_HIGHORDER_TOPOLOGY_API GenerateBezierCylinder<Vec3fTypes>;
#endif //SOFA_DOUBLE


} // namespace constraint

} // namespace component

} // namespace sofa

