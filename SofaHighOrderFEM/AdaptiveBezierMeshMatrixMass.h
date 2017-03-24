#ifndef SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERMESHMATRIXMASS_H
#define SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERMESHMATRIXMASS_H


#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif
#include "HighOrderMeshMatrixMass.h"


namespace sofa
{
namespace component
{
namespace topology 
{
	class AdaptiveHighOrderTetrahedronSetTopologyContainer;
}
namespace mass
{



using namespace sofa::defaulttype;
using namespace sofa::component::topology;

template <class DataTypes, class TMassType>
class AdaptiveBezierHighOrderMeshMatrixMass : public HighOrderMeshMatrixMass<DataTypes,TMassType>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(AdaptiveBezierHighOrderMeshMatrixMass,DataTypes,TMassType), SOFA_TEMPLATE2(HighOrderMeshMatrixMass,DataTypes,TMassType));

    typedef HighOrderMeshMatrixMass<DataTypes,TMassType> Inherited;
    typedef typename DataTypes::VecCoord                    VecCoord;
    typedef typename DataTypes::VecDeriv                    VecDeriv;
    typedef typename DataTypes::Coord                       Coord;
    typedef typename DataTypes::Deriv                       Deriv;
    typedef typename DataTypes::Real                        Real;
    typedef core::objectmodel::Data<VecCoord>               DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv>               DataVecDeriv;
    typedef TMassType                                       MassType;
    typedef helper::vector<MassType>                        MassVector;
    typedef helper::vector<MassVector>                      MassVectorVector;
	typedef sofa::core::topology::BaseMeshTopology::TetrahedraAroundTriangle TetrahedraAroundTriangle;
	typedef sofa::core::topology::BaseMeshTopology::TetrahedraAroundEdge TetrahedraAroundEdge;
	typedef sofa::core::topology::BaseMeshTopology::TetrahedraAroundVertex TetrahedraAroundVertex;

    // In case of non 3D template
    typedef defaulttype::Vec<3,Real> Vec3;

    /// assumes the geometry object type is 3D
    typedef defaulttype::StdVectorTypes< Vec3, Vec3, Real > GeometricalTypes ;

protected:  
    sofa::component::topology::AdaptiveHighOrderTetrahedronSetTopologyContainer *adaptiveContainer;

    /* ---------- Specific data for Bezier Elements ------*/
    /// use this data structure to store mass for Bezier tetrahedra. 
    //// The size of the vector is nbControlPoints*(nbControlPoints+1)/2 where nbControlPoints=(degree+1)*(degree+2)*(degree+3)/2
	topology::TetrahedronData<helper::vector<MassVector> > reducedTetrahedronMassInfo;
	AdaptiveBezierHighOrderMeshMatrixMass();

	virtual ~AdaptiveBezierHighOrderMeshMatrixMass();
 
public:


	virtual void clear();

	virtual void handleTopologyChange(core::topology::Topology *topo);
	virtual void checkConsistencyOfLumpedMass();
    virtual void init();

   

	// returns the mass vector for a given index of a Bezier tetrahedron
	virtual const  MassVector &getBezierTetrahedronMassVector(const size_t i) const;



};



} // namespace mass

} // namespace component

} // namespace sofa

#endif
