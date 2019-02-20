
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
        .add< BezierTetrahedronSetGeometryAlgorithms<Vec3Types> >(true) // default template
        .add< BezierTetrahedronSetGeometryAlgorithms<Vec2Types> >()
        .add< BezierTetrahedronSetGeometryAlgorithms<Vec1Types> >();


template class SOFA_HIGHORDER_TOPOLOGY_API BezierTetrahedronSetGeometryAlgorithms<Vec3Types>;
template class SOFA_HIGHORDER_TOPOLOGY_API BezierTetrahedronSetGeometryAlgorithms<Vec2Types>;
template class SOFA_HIGHORDER_TOPOLOGY_API BezierTetrahedronSetGeometryAlgorithms<Vec1Types>;




} // namespace topology

} // namespace component

} // namespace sofa

