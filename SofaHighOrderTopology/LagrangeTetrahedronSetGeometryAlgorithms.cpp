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
#ifdef SOFA_FLOAT
        .add< LagrangeTetrahedronSetGeometryAlgorithms<Vec3fTypes> >(true) // default template
#else
        .add< LagrangeTetrahedronSetGeometryAlgorithms<Vec3dTypes> >(true) // default template
#ifndef SOFA_DOUBLE
        .add< LagrangeTetrahedronSetGeometryAlgorithms<Vec3fTypes> >() // default template
#endif
#endif
#ifndef SOFA_FLOAT
        .add< LagrangeTetrahedronSetGeometryAlgorithms<Vec2dTypes> >()
        .add< LagrangeTetrahedronSetGeometryAlgorithms<Vec1dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< LagrangeTetrahedronSetGeometryAlgorithms<Vec2fTypes> >()
        .add< LagrangeTetrahedronSetGeometryAlgorithms<Vec1fTypes> >()
#endif
        ;

#ifndef SOFA_FLOAT
template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTetrahedronSetGeometryAlgorithms<Vec3dTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTetrahedronSetGeometryAlgorithms<Vec2dTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTetrahedronSetGeometryAlgorithms<Vec1dTypes>;
#endif

#ifndef SOFA_DOUBLE
template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTetrahedronSetGeometryAlgorithms<Vec3fTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTetrahedronSetGeometryAlgorithms<Vec2fTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTetrahedronSetGeometryAlgorithms<Vec1fTypes>;
#endif

} // namespace topology

} // namespace component

} // namespace sofa

