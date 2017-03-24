
#include "AdaptiveBezierTetrahedronSetTopologyAlgorithms.inl"
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace topology
{

using namespace sofa::defaulttype;
SOFA_DECL_CLASS(AdaptiveBezierTetrahedronSetTopologyAlgorithms)
int AdaptiveBezierTetrahedronSetTopologyAlgorithmsClass = core::RegisterObject("Adaptive Bezier Tetrahedron set topology algorithms")
#ifdef SOFA_FLOAT
        .add< AdaptiveBezierTetrahedronSetTopologyAlgorithms<Vec3fTypes> >(true) // default template
#else
        .add< AdaptiveBezierTetrahedronSetTopologyAlgorithms<Vec3dTypes> >(true) // default template
#ifndef SOFA_DOUBLE
        .add< AdaptiveBezierTetrahedronSetTopologyAlgorithms<Vec3fTypes> >() // default template
#endif
#endif
        ;

#ifndef SOFA_FLOAT
template class SOFA_HIGHORDER_FEM_API AdaptiveBezierTetrahedronSetTopologyAlgorithms<Vec3dTypes>;

#endif

#ifndef SOFA_DOUBLE
template class SOFA_HIGHORDER_FEM_API AdaptiveBezierTetrahedronSetTopologyAlgorithms<Vec3fTypes>;

#endif

} // namespace topology

} // namespace component

} // namespace sofa

