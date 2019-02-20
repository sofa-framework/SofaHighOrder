#define SOFA_HIGHORDERTOPOLOGY_LAGRANGETRIANGLESETGEOMETRYALGORITHMS_CPP
#include "LagrangeTriangleSetGeometryAlgorithms.inl"
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



SOFA_DECL_CLASS(LagrangeTriangleSetGeometryAlgorithms)
int LagrangeTriangleSetGeometryAlgorithmsClass = core::RegisterObject("Lagrange Triangle set geometry algorithms")
        .add< LagrangeTriangleSetGeometryAlgorithms<Vec3Types> >(true) // default template
        .add< LagrangeTriangleSetGeometryAlgorithms<Vec2Types> >()
        .add< LagrangeTriangleSetGeometryAlgorithms<Vec1Types> >()
        ;


template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTriangleSetGeometryAlgorithms<Vec3Types>;
template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTriangleSetGeometryAlgorithms<Vec2Types>;
template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTriangleSetGeometryAlgorithms<Vec1Types>;


} // namespace topology

} // namespace component

} // namespace sofa

