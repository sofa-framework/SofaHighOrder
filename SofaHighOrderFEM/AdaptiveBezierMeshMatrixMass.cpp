#define SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERMESHMATRIXMASS_CPP
#include "AdaptiveBezierMeshMatrixMass.inl"
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/gl/Axis.h>
#include "initHighOrderFEM.h"

namespace sofa
{

namespace component
{

namespace mass
{

using namespace sofa::defaulttype;




SOFA_DECL_CLASS(AdaptiveBezierHighOrderMeshMatrixMass)

// Register in the Factory
int AdaptiveBezierHighOrderMeshMatrixMassClass = core::RegisterObject("Define a mass matrix for an adaptive Bezier mesh")
#ifndef SOFA_FLOAT
        .add< AdaptiveBezierHighOrderMeshMatrixMass<Vec3dTypes,double> >()

#endif
#ifndef SOFA_DOUBLE
        .add< AdaptiveBezierHighOrderMeshMatrixMass<Vec3fTypes,float> >()

#endif
        ;

#ifndef SOFA_FLOAT
template class SOFA_HIGHORDER_FEM_API AdaptiveBezierHighOrderMeshMatrixMass<Vec3dTypes,double>;

#endif
#ifndef SOFA_DOUBLE
template class SOFA_HIGHORDER_FEM_API AdaptiveBezierHighOrderMeshMatrixMass<Vec3fTypes,float>;

#endif


} // namespace mass

} // namespace component

} // namespace sofa

