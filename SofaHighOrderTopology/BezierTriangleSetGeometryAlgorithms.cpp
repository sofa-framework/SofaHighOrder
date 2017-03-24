
#define SOFA_HIGHORDERTOPOLOGY_BEZIERTRIANGLESETGEOMETRYALGORITHMS_CPP
#include "BezierTriangleSetGeometryAlgorithms.inl"
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
SOFA_DECL_CLASS(BezierTriangleSetGeometryAlgorithms)
int BezierTriangleSetGeometryAlgorithmsClass = core::RegisterObject("Bezier Triangle set geometry algorithms")
#ifdef SOFA_FLOAT
        .add< BezierTriangleSetGeometryAlgorithms<Vec3fTypes> >(true) // default template
#else
        .add< BezierTriangleSetGeometryAlgorithms<Vec3dTypes> >(true) // default template
#ifndef SOFA_DOUBLE
        .add< BezierTriangleSetGeometryAlgorithms<Vec3fTypes> >() // default template
#endif
#endif
#ifndef SOFA_FLOAT
        .add< BezierTriangleSetGeometryAlgorithms<Vec2dTypes> >()
        .add< BezierTriangleSetGeometryAlgorithms<Vec1dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< BezierTriangleSetGeometryAlgorithms<Vec2fTypes> >()
        .add< BezierTriangleSetGeometryAlgorithms<Vec1fTypes> >()
#endif
        ;

#ifndef SOFA_FLOAT
template class SOFA_HIGHORDER_TOPOLOGY_API BezierTriangleSetGeometryAlgorithms<Vec3dTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API BezierTriangleSetGeometryAlgorithms<Vec2dTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API BezierTriangleSetGeometryAlgorithms<Vec1dTypes>;
#endif

#ifndef SOFA_DOUBLE
template class SOFA_HIGHORDER_TOPOLOGY_API BezierTriangleSetGeometryAlgorithms<Vec3fTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API BezierTriangleSetGeometryAlgorithms<Vec2fTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API BezierTriangleSetGeometryAlgorithms<Vec1fTypes>;
#endif

} // namespace topology

} // namespace component

} // namespace sofa

