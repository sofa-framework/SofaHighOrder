
#define SOFA_HIGHORDERTOPOLOGY_BEZIER2MESHMECHANICALMAPPING_CPP
#include "Bezier2MeshMechanicalMapping.inl"


#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace mapping
{

using namespace sofa::defaulttype;


SOFA_DECL_CLASS(Bezier2MeshMechanicalMapping)

int Bezier2MeshMechanicalMappingClass = core::RegisterObject("Mechanical mapping between a Bezier triangle or Bezier tetrahedra with a tesselated triangle mesh or tesselated tetrahedron mesh")
        .add< Bezier2MeshMechanicalMapping<Vec3Types, Vec3Types> >(true) // default template
        .add< Bezier2MeshMechanicalMapping<Vec2Types,Vec2Types> >()
        .add< Bezier2MeshMechanicalMapping<Vec1Types,Vec1Types> >();


template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< Vec3Types, Vec3Types >;
template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< Vec2Types, Vec2Types >;
template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< Vec1Types, Vec1Types >;
template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< Vec3Types, ExtVec3Types >;



} // namespace mapping

} // namespace component

} // namespace sofa
