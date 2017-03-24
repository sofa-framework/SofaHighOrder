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
#ifdef SOFA_FLOAT
        .add< LagrangeTriangleSetGeometryAlgorithms<Vec3fTypes> >(true) // default template
#else
        .add< LagrangeTriangleSetGeometryAlgorithms<Vec3dTypes> >(true) // default template
#ifndef SOFA_DOUBLE
        .add< LagrangeTriangleSetGeometryAlgorithms<Vec3fTypes> >() // default template
#endif
#endif
#ifndef SOFA_FLOAT
        .add< LagrangeTriangleSetGeometryAlgorithms<Vec2dTypes> >()
        .add< LagrangeTriangleSetGeometryAlgorithms<Vec1dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        .add< LagrangeTriangleSetGeometryAlgorithms<Vec2fTypes> >()
        .add< LagrangeTriangleSetGeometryAlgorithms<Vec1fTypes> >()
#endif
        ;

#ifndef SOFA_FLOAT
template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTriangleSetGeometryAlgorithms<Vec3dTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTriangleSetGeometryAlgorithms<Vec2dTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTriangleSetGeometryAlgorithms<Vec1dTypes>;
#endif

#ifndef SOFA_DOUBLE
template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTriangleSetGeometryAlgorithms<Vec3fTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTriangleSetGeometryAlgorithms<Vec2fTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTriangleSetGeometryAlgorithms<Vec1fTypes>;
#endif

} // namespace topology

} // namespace component

} // namespace sofa

