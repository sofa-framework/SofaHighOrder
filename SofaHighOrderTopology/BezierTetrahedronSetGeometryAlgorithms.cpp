
#define SOFA_HIGHORDERTOPOLOGY_BEZIERTETRAHEDRONSETGEOMETRYALGORITHMS_CPP
#include "BezierTetrahedronSetGeometryAlgorithms.inl"
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/core/ObjectFactory.h>
#include "initHighOrderTopology.h"
namespace sofa
{

namespace component
{

namespace topology
{
using namespace sofa::defaulttype;
SOFA_DECL_CLASS(BezierTetrahedronSetGeometryAlgorithms)
int BezierTetrahedronSetGeometryAlgorithmsClass = core::RegisterObject("Bezier Tetrahedron set geometry algorithms")
#ifdef SOFA_FLOAT
        .add< BezierTetrahedronSetGeometryAlgorithms<Vec3fTypes> >(true) // default template
#else
        .add< BezierTetrahedronSetGeometryAlgorithms<Vec3dTypes> >(true) // default template
#ifndef SOFA_DOUBLE
        .add< BezierTetrahedronSetGeometryAlgorithms<Vec3fTypes> >() // default template
#endif
#endif
#ifndef SOFA_FLOAT
        .add< BezierTetrahedronSetGeometryAlgorithms<Vec2dTypes> >()
        .add< BezierTetrahedronSetGeometryAlgorithms<Vec1dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< BezierTetrahedronSetGeometryAlgorithms<Vec2fTypes> >()
        .add< BezierTetrahedronSetGeometryAlgorithms<Vec1fTypes> >()
#endif
        ;

#ifndef SOFA_FLOAT
template class SOFA_HIGHORDER_TOPOLOGY_API BezierTetrahedronSetGeometryAlgorithms<Vec3dTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API BezierTetrahedronSetGeometryAlgorithms<Vec2dTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API BezierTetrahedronSetGeometryAlgorithms<Vec1dTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API HighOrderTetrahedronSetGeometryAlgorithms<Vec3dTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API HighOrderTetrahedronSetGeometryAlgorithms<Vec2dTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API HighOrderTetrahedronSetGeometryAlgorithms<Vec1dTypes>;
#endif

#ifndef SOFA_DOUBLE
template class SOFA_HIGHORDER_TOPOLOGY_API BezierTetrahedronSetGeometryAlgorithms<Vec3fTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API BezierTetrahedronSetGeometryAlgorithms<Vec2fTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API BezierTetrahedronSetGeometryAlgorithms<Vec1fTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API HighOrderTetrahedronSetGeometryAlgorithms<Vec3fTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API HighOrderTetrahedronSetGeometryAlgorithms<Vec2fTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API HighOrderTetrahedronSetGeometryAlgorithms<Vec1fTypes>;

#endif

} // namespace topology

} // namespace component

} // namespace sofa

