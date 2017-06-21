
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
#ifdef SOFA_FLOAT
        .add< Bezier2MeshMechanicalMapping<Vec3fTypes, Vec3fTypes> >(true) // default template
#else
        .add< Bezier2MeshMechanicalMapping<Vec3dTypes, Vec3dTypes> >(true) // default template
#ifndef SOFA_DOUBLE
        .add< Bezier2MeshMechanicalMapping<Vec3fTypes,Vec3fTypes> >() // default template
#endif
#endif
#ifndef SOFA_FLOAT
        .add< Bezier2MeshMechanicalMapping<Vec2dTypes,Vec2dTypes> >()
        .add< Bezier2MeshMechanicalMapping<Vec1dTypes,Vec1dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< Bezier2MeshMechanicalMapping<Vec2fTypes,Vec2fTypes> >()
        .add< Bezier2MeshMechanicalMapping<Vec1fTypes,Vec1fTypes> >()
#endif
        ;


#ifndef SOFA_FLOAT
template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< Vec3dTypes, Vec3dTypes >;
template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< Vec2dTypes, Vec2dTypes >;
template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< Vec1dTypes, Vec1dTypes >;
template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< Vec3dTypes, ExtVec3dTypes >;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< Vec3fTypes, Vec3fTypes >;
template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< Vec2fTypes, Vec2fTypes >;
template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< Vec3fTypes, ExtVec3fTypes >;
template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< Vec1fTypes, Vec1fTypes >;
#endif

#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< Vec3dTypes, Vec3fTypes >;
template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< Vec3fTypes, Vec3dTypes >;
template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< Vec2fTypes, Vec2dTypes >;
template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< Vec2dTypes, Vec2fTypes >;
template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< Vec1dTypes, Vec1fTypes >;
template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< Vec1fTypes, Vec1dTypes >;
template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< Vec3fTypes, ExtVec3dTypes >;
template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< Vec3dTypes, ExtVec3fTypes >;
#endif
#endif

} // namespace mapping

} // namespace component

} // namespace sofa
