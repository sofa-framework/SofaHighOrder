#define SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERTETRA2BEZIERTRIANGLEMECHANICALMAPPING_CPP
#include "AdaptiveBezierTetra2BezierTriangleMechanicalMapping.inl"


#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace mapping
{

using namespace sofa::defaulttype;


SOFA_DECL_CLASS(AdaptiveBezierTetra2BezierTriangleMechanicalMapping)

int AdaptiveBezierTetra2BezierTriangleMechanicalMappingClass = core::RegisterObject("Mechanical mapping between an Adaptive Bezier tetrahedralization and a Bezier Triangulation")
#ifndef SOFA_FLOAT
        .add< AdaptiveBezierTetra2BezierTriangleMechanicalMapping< Vec3dTypes, Vec3dTypes > >()
        .add< AdaptiveBezierTetra2BezierTriangleMechanicalMapping< Vec3dTypes, ExtVec3dTypes > >()
#endif
#ifndef SOFA_DOUBLE
        .add< AdaptiveBezierTetra2BezierTriangleMechanicalMapping< Vec3fTypes, Vec3fTypes > >()
        .add< AdaptiveBezierTetra2BezierTriangleMechanicalMapping< Vec3fTypes, ExtVec3fTypes > >()
#endif

#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
        .add< AdaptiveBezierTetra2BezierTriangleMechanicalMapping< Vec3fTypes, Vec3dTypes > >()
        .add< AdaptiveBezierTetra2BezierTriangleMechanicalMapping< Vec3dTypes, Vec3fTypes > >()
		.add< AdaptiveBezierTetra2BezierTriangleMechanicalMapping< Vec3dTypes, ExtVec3fTypes > >()
		.add< AdaptiveBezierTetra2BezierTriangleMechanicalMapping< Vec3fTypes, ExtVec3dTypes > >()
#endif
#endif
//.addAlias("SimpleTesselatedTetraMechanicalMapping")
        ;


#ifndef SOFA_FLOAT

template class SOFA_HIGHORDER_FEM_API AdaptiveBezierTetra2BezierTriangleMechanicalMapping< Vec3dTypes, Vec3dTypes >;
//template class SOFA_HIGHORDER_FEM_API AdaptiveBezierTetra2BezierTriangleMechanicalMapping< Vec3dTypes, Vec3dTypes >;
template class SOFA_HIGHORDER_FEM_API AdaptiveBezierTetra2BezierTriangleMechanicalMapping< Vec3dTypes, ExtVec3dTypes >;
#endif
#ifndef SOFA_DOUBLE

template class SOFA_HIGHORDER_FEM_API AdaptiveBezierTetra2BezierTriangleMechanicalMapping< Vec3fTypes, Vec3fTypes >;
template class SOFA_HIGHORDER_FEM_API AdaptiveBezierTetra2BezierTriangleMechanicalMapping< Vec3fTypes, ExtVec3fTypes >;
#endif

#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_HIGHORDER_FEM_API AdaptiveBezierTetra2BezierTriangleMechanicalMapping< Vec3dTypes, Vec3fTypes >;
template class SOFA_HIGHORDER_FEM_API AdaptiveBezierTetra2BezierTriangleMechanicalMapping< Vec3fTypes, Vec3dTypes >;
template class SOFA_HIGHORDER_FEM_API AdaptiveBezierTetra2BezierTriangleMechanicalMapping< Vec3fTypes, ExtVec3dTypes >;
template class SOFA_HIGHORDER_FEM_API AdaptiveBezierTetra2BezierTriangleMechanicalMapping< Vec3dTypes, ExtVec3fTypes >;
#endif
#endif

} // namespace mapping
namespace topology
{
template class SOFA_HIGHORDER_FEM_API PointData< sofa::helper::vector<sofa::defaulttype::Vec3dTypes::Coord >  >;
template class SOFA_HIGHORDER_FEM_API PointData< sofa::helper::vector<sofa::defaulttype::Vec3fTypes::Coord >  >;
}
} // namespace component

} // namespace sofa
