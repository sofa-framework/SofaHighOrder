#define SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERTETRAHEDRONCONTROLLER_CPP

#include "AdaptiveBezierTetrahedronController.inl"
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace controller
{

SOFA_DECL_CLASS(AdaptiveBezierTetrahedronController)

using namespace sofa::defaulttype;

// Register in the Factory
int AdaptiveBezierTetrahedronControllerClass = core::RegisterObject("Control the degree of Bezier tetrahedra on AdaptiveBezierTetrahedronContainer")
#ifndef SOFA_FLOAT
        .add< AdaptiveBezierTetrahedronController<Vec3dTypes> >(true)
#endif
#ifndef SOFA_DOUBLE
        .add< AdaptiveBezierTetrahedronController<Vec3fTypes> >()

#endif
        ;



#ifndef SOFA_FLOAT
template class SOFA_HIGHORDER_FEM_API AdaptiveBezierTetrahedronController<Vec3dTypes>;

#endif
#ifndef SOFA_DOUBLE
template class  SOFA_HIGHORDER_FEM_API AdaptiveBezierTetrahedronController<Vec3fTypes>;

#endif


} // namespace misc

} // namespace component

} // namespace sofa
