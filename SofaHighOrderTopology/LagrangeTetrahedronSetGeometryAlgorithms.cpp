#define SOFA_HIGHORDERTOPOLOGY_LAGRANGETETRAHEDRONSETGEOMETRYALGORITHMS_CPP
#include "LagrangeTetrahedronSetGeometryAlgorithms.inl"
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
SOFA_DECL_CLASS(LagrangeTetrahedronSetGeometryAlgorithms)
int LagrangeTetrahedronSetGeometryAlgorithmsClass = core::RegisterObject("Lagrange Tetrahedron set geometry algorithms")
        .add< LagrangeTetrahedronSetGeometryAlgorithms<Vec3Types> >(true) // default template
        .add< LagrangeTetrahedronSetGeometryAlgorithms<Vec2Types> >()
        ;

template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTetrahedronSetGeometryAlgorithms<Vec3Types>;
template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTetrahedronSetGeometryAlgorithms<Vec2Types>;
template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTetrahedronSetGeometryAlgorithms<Vec1Types>;

} // namespace topology

} // namespace component

} // namespace sofa

