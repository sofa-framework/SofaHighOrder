
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



template <>
Vector6 HighOrderMeshMatrixMass<Vec3Types, double>::getMomentum ( const core::MechanicalParams*, const DataVecCoord& vx, const DataVecDeriv& vv ) const
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








SOFA_DECL_CLASS(HighOrderMeshMatrixMass)

// Register in the Factory
int HighOrderMeshMatrixMassClass = core::RegisterObject("Define a mass matrix for high order tetrahedron elements")
        .add< HighOrderMeshMatrixMass<Vec3Types, Vec3Types::Real> >()
	    .add< HighOrderMeshMatrixMass<Vec2Types, Vec2Types::Real> >()
	    .add< HighOrderMeshMatrixMass<Vec1Types, Vec1Types::Real> >()
        ;


template class SOFA_HIGHORDER_FEM_API HighOrderMeshMatrixMass<Vec3Types, Vec3Types::Real>;
template class SOFA_HIGHORDER_FEM_API HighOrderMeshMatrixMass<Vec2Types, Vec2Types::Real>;
template class SOFA_HIGHORDER_FEM_API HighOrderMeshMatrixMass<Vec1Types, Vec1Types::Real>;


} // namespace mass

} // namespace component

} // namespace sofa

