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
        .add< GenerateBezierCylinder<Vec3Types> >()
        ;

template class SOFA_HIGHORDER_TOPOLOGY_API GenerateBezierCylinder<Vec3Types>;


} // namespace constraint

} // namespace component

} // namespace sofa

