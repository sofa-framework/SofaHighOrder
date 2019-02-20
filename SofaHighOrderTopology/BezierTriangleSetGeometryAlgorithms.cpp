
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
        .add< BezierTriangleSetGeometryAlgorithms<Vec3Types> >(true) // default template
        .add< BezierTriangleSetGeometryAlgorithms<Vec2Types> >()
        .add< BezierTriangleSetGeometryAlgorithms<Vec1Types> >();

template class SOFA_HIGHORDER_TOPOLOGY_API BezierTriangleSetGeometryAlgorithms<Vec3Types>;
template class SOFA_HIGHORDER_TOPOLOGY_API BezierTriangleSetGeometryAlgorithms<Vec2Types>;
template class SOFA_HIGHORDER_TOPOLOGY_API BezierTriangleSetGeometryAlgorithms<Vec1Types>;


} // namespace topology

} // namespace component

} // namespace sofa

