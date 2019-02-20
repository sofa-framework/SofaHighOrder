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
        .add< GenerateBezierSphere<Vec3Types> >()
        ;

template class SOFA_HIGHORDER_TOPOLOGY_API GenerateBezierSphere<Vec3Types>;


} // namespace constraint

} // namespace component

} // namespace sofa
