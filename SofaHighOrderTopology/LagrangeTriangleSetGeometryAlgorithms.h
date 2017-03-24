#ifndef SOFA_HIGHORDERTOPOLOGY_LAGRANGETRIANGLESETGEOMETRYALGORITHMS_H
#define SOFA_HIGHORDERTOPOLOGY_LAGRANGETRIANGLESETGEOMETRYALGORITHMS_H


#include <SofaBaseTopology/TriangleSetGeometryAlgorithms.h>
#include "HighOrderTriangleSetTopologyContainer.h"
#include "HighOrderTriangleSetGeometryAlgorithms.h"

namespace sofa
{

namespace component
{

namespace topology
{
using core::topology::BaseMeshTopology;


/**
* A class that provides geometry information on a high order Lagrange Tetrahedral Mesh
*/
template < class DataTypes >
class LagrangeTriangleSetGeometryAlgorithms : public HighOrderTriangleSetGeometryAlgorithms<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(LagrangeTriangleSetGeometryAlgorithms,DataTypes),SOFA_TEMPLATE(HighOrderTriangleSetGeometryAlgorithms,DataTypes));
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef sofa::defaulttype::Vec<3,Real> Vec3;
    typedef sofa::defaulttype::Mat<3,3,Real> Mat33;

	typedef BaseMeshTopology::TetraID TetraID;
	typedef BaseMeshTopology::Tetra Tetra;
	typedef BaseMeshTopology::SeqTetrahedra SeqTetrahedra;
	typedef BaseMeshTopology::SeqEdges SeqEdges;
	typedef BaseMeshTopology::TetrahedraAroundVertex TetrahedraAroundVertex;
	typedef BaseMeshTopology::TetrahedraAroundEdge TetrahedraAroundEdge;
	typedef BaseMeshTopology::TetrahedraAroundTriangle TetrahedraAroundTriangle;
	typedef BaseMeshTopology::EdgesInTriangle EdgesInTriangle;
	typedef HighOrderTriangleSetTopologyContainer::VecPointID VecPointID;

	typedef Tetra Triangle;
protected:

	/// constructor 
	LagrangeTriangleSetGeometryAlgorithms();
    virtual ~LagrangeTriangleSetGeometryAlgorithms() {}
public:
	virtual void init();

	/// computes the shape function 
	virtual Real computeShapeFunction(const TriangleIndexVector tbi, const Vec3 barycentricCoordinate);
	/// computes the shape function gradient
    virtual Vec3 computeShapeFunctionDerivatives(const TriangleIndexVector tbi, const Vec3 barycentricCoordinate);
    /// computes the shape function hessian
    virtual Mat33 computeShapeFunctionHessian(const TriangleIndexVector tbi, const Vec3 barycentricCoordinate);

	/// returns the mass coefficient used to compute the mass matrix when a triangle is affine 
	virtual Real getAffineMassCoefficient(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2);
	/// returns the mass coefficient used to compute the exact mass matrix for  a regular high order element  
	virtual Real getRegularMassCoefficient(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2,const size_t indr[2]);
	/// returns the stiffness coefficient used to compute the stiffness matrix when a triangle is affine 
	virtual Mat33 getAffineStiffnessCoefficientMatrix(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2);
private:
	helper::vector< helper::vector<Real> > stirlingNumberArray;
	void checkCoefficients();	
	long int sterlingVector(const TriangleIndexVector tbiIn,const TriangleIndexVector tbiIk) const;
	/// compute  the mass coefficient used for  the exact mass matrix for  a regular high order element  
	 Real computeRegularMassCoefficient(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2,const size_t indr[2]);
	/// compute the mass coefficient for  the mass matrix when a triangle is affine 
	 Real computeAffineMassCoefficient(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2);
	/// compute the stiffness coefficient for the stiffness matrix when a triangle is affine 
	Mat33 computeAffineStiffnessCoefficientMatrix(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2);

};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_HIGHORDERTOPOLOGY_LAGRANGETRIANGLESETGEOMETRYALGORITHMS_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTriangleSetGeometryAlgorithms<defaulttype::Vec3dTypes>;
extern template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTriangleSetGeometryAlgorithms<defaulttype::Vec2dTypes>;
extern template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTriangleSetGeometryAlgorithms<defaulttype::Vec1dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTriangleSetGeometryAlgorithms<defaulttype::Vec3fTypes>;
extern template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTriangleSetGeometryAlgorithms<defaulttype::Vec2fTypes>;
extern template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTriangleSetGeometryAlgorithms<defaulttype::Vec1fTypes>;
#endif
#endif

} // namespace topology

} // namespace component

} // namespace sofa

#endif
