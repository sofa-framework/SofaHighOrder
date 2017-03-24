#define SOFA_HIGHORDERTOPOLOGY_GENERATEBEZIERSPHERE_CPP
#include "GenerateBezierSphere.inl"
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/Vec3Types.h>

namespace sofa
{

namespace component
{

namespace engine
{
using namespace sofa::defaulttype;

SOFA_DECL_CLASS(GenerateBezierSphere)

int GenerateBezierSphereClass = core::RegisterObject("Generate a sphereical (Bezier) Tetrahedral and Triangular Mesh")
#ifndef SOFA_FLOAT
        .add< GenerateBezierSphere<Vec3dTypes> >()
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
        .add< GenerateBezierSphere<Vec3fTypes> >()
#endif //SOFA_DOUBLE
        ;


#ifndef SOFA_FLOAT
template class SOFA_HIGHORDER_TOPOLOGY_API GenerateBezierSphere<Vec3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_HIGHORDER_TOPOLOGY_API GenerateBezierSphere<Vec3fTypes>;
#endif //SOFA_DOUBLE


} // namespace constraint

} // namespace component

} // namespace sofa
