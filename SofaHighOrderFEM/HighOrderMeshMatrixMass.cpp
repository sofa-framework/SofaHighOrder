
#define SOFA_HIGHORDERDFEM_MESHMATRIXMASS_CPP
#include "HighOrderMeshMatrixMass.inl"
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/gl/Axis.h>

namespace sofa
{

namespace component
{

namespace mass
{

using namespace sofa::defaulttype;


#ifndef SOFA_FLOAT

template <>
Vector6 HighOrderMeshMatrixMass<Vec3dTypes, double>::getMomentum ( const core::MechanicalParams*, const DataVecCoord& vx, const DataVecDeriv& vv ) const
{
    const MassVector &vertexMass= vertexMassInfo.getValue();
    const MassVector &edgeMass= edgeMassInfo.getValue();

    helper::ReadAccessor< DataVecCoord > x = vx;
    helper::ReadAccessor< DataVecDeriv > v = vv;

    Vector6 momentum;
    for( unsigned int i=0 ; i<v.size() ; i++ )
    {
        Deriv linearMomentum = v[i] * vertexMass[i];
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[j] += linearMomentum[j];
        Deriv angularMomentum = cross( x[i], linearMomentum );
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[3+j] += angularMomentum[j];
    }

    for( int i=0 ; i<_topology->getNbEdges() ; ++i )
    {
        unsigned v0 = _topology->getEdge(i)[0];
        unsigned v1 = _topology->getEdge(i)[1];

        // is it correct to share the edge mass between the 2 vertices?
        double m = edgeMass[i] * 0.5;

        Deriv linearMomentum = v[v0] * m;
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[j] += linearMomentum[j];
        Deriv angularMomentum = cross( x[v0], linearMomentum );
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[3+j] += angularMomentum[j];

        linearMomentum = v[v1] * m;
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[j] += linearMomentum[j];
        angularMomentum = cross( x[v1], linearMomentum );
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[3+j] += angularMomentum[j];
    }

    return momentum;
}

#endif
#ifndef SOFA_DOUBLE

template <>
Vector6 HighOrderMeshMatrixMass<Vec3fTypes, float>::getMomentum ( const core::MechanicalParams*, const DataVecCoord& vx, const DataVecDeriv& vv ) const
{
    const MassVector &vertexMass= vertexMassInfo.getValue();
    const MassVector &edgeMass= edgeMassInfo.getValue();

    helper::ReadAccessor< DataVecCoord > x = vx;
    helper::ReadAccessor< DataVecDeriv > v = vv;

    Vector6 momentum;
    for( unsigned int i=0 ; i<v.size() ; i++ )
    {
        Deriv linearMomentum = v[i] * vertexMass[i];
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[j] += linearMomentum[j];
        Deriv angularMomentum = cross( x[i], linearMomentum );
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[3+j] += angularMomentum[j];
    }

    for( int i=0 ; i<_topology->getNbEdges() ; ++i )
    {
        unsigned v0 = _topology->getEdge(i)[0];
        unsigned v1 = _topology->getEdge(i)[1];

        // is it correct to share the edge mass between the 2 vertices?
        float m = edgeMass[i] * 0.5f;

        Deriv linearMomentum = v[v0] * m;
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[j] += linearMomentum[j];
        Deriv angularMomentum = cross( x[v0], linearMomentum );
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[3+j] += angularMomentum[j];

        linearMomentum = v[v1] * m;
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[j] += linearMomentum[j];
        angularMomentum = cross( x[v1], linearMomentum );
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[3+j] += angularMomentum[j];
    }

    return momentum;
}


#endif







SOFA_DECL_CLASS(HighOrderMeshMatrixMass)

// Register in the Factory
int HighOrderMeshMatrixMassClass = core::RegisterObject("Define a mass matrix for high order tetrahedron elements")
#ifndef SOFA_FLOAT
        .add< HighOrderMeshMatrixMass<Vec3dTypes,double> >()
	    .add< HighOrderMeshMatrixMass<Vec2dTypes,double> >()
	    .add< HighOrderMeshMatrixMass<Vec1dTypes,double> >()

#endif
#ifndef SOFA_DOUBLE
        .add< HighOrderMeshMatrixMass<Vec3fTypes,float> >()
   	    .add< HighOrderMeshMatrixMass<Vec2fTypes,float> >()
   	    .add< HighOrderMeshMatrixMass<Vec1fTypes,float> >()

#endif
        ;

#ifndef SOFA_FLOAT
template class SOFA_HIGHORDER_FEM_API HighOrderMeshMatrixMass<Vec3dTypes,double>;
template class SOFA_HIGHORDER_FEM_API HighOrderMeshMatrixMass<Vec2dTypes,double>;
template class SOFA_HIGHORDER_FEM_API HighOrderMeshMatrixMass<Vec1dTypes,double>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_HIGHORDER_FEM_API HighOrderMeshMatrixMass<Vec3fTypes,float>;
template class SOFA_HIGHORDER_FEM_API HighOrderMeshMatrixMass<Vec2fTypes,float>;
template class SOFA_HIGHORDER_FEM_API HighOrderMeshMatrixMass<Vec1fTypes,float>;
#endif


} // namespace mass

} // namespace component

} // namespace sofa

