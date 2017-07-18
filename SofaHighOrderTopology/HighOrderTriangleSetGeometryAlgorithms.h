#ifndef SOFA_HIGHORDERTOPOLOGY_HIGHORDERTRIANGLESETGEOMETRYALGORITHMS_H
#define SOFA_HIGHORDERTOPOLOGY_HIGHORDERTRIANGLESETGEOMETRYALGORITHMS_H


#include <SofaBaseTopology/TriangleSetGeometryAlgorithms.h>
#include "HighOrderTriangleSetTopologyContainer.h"


namespace sofa
{

namespace component
{

namespace topology
{
using core::topology::BaseMeshTopology;


template< typename Real, int N> typename NumericalIntegrationDescriptor<Real, N>::QuadraturePointArray TriangleConicalRule(const size_t d);
/**
* A class that provides geometry information on a High Order Triangle mesh 
*/
template < class DataTypes >
class HighOrderTriangleSetGeometryAlgorithms : public TriangleSetGeometryAlgorithms<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(HighOrderTriangleSetGeometryAlgorithms,DataTypes),SOFA_TEMPLATE(TriangleSetGeometryAlgorithms,DataTypes));
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef sofa::defaulttype::Vec<3,Real> Vec3;
    typedef sofa::defaulttype::Mat<3,3,Real> Mat33;

	typedef BaseMeshTopology::Edge Edge;
	typedef BaseMeshTopology::PointID PointID;
	typedef BaseMeshTopology::TetraID TetraID;
	typedef BaseMeshTopology::Triangle Triangle;
	typedef BaseMeshTopology::SeqTetrahedra SeqTetrahedra;
	typedef BaseMeshTopology::SeqEdges SeqEdges;
	typedef BaseMeshTopology::TetrahedraAroundVertex TetrahedraAroundVertex;
	typedef BaseMeshTopology::TetrahedraAroundEdge TetrahedraAroundEdge;
	typedef BaseMeshTopology::TetrahedraAroundTriangle TetrahedraAroundTriangle;
	typedef BaseMeshTopology::EdgesInTriangle EdgesInTriangle;
	typedef HighOrderTriangleSetTopologyContainer::VecPointID VecPointID;

protected:
   
	/// container	
	HighOrderTriangleSetTopologyContainer *container; 
	/// degree of the polynomial
	HighOrderDegreeType degree; 
	// array of Tetrahedral Bezier indices
	sofa::helper::vector<TriangleIndexVector> tbiArray;
	/// the list of edges of the Bezier Triangle used in the draw function
    std::set<std::pair<Edge,size_t> > bezierTriangleEdgeSet;
    // if the new cubature tables have been added
    bool initializedNewCubatureTables;
    // add new cubature rules
    void defineNewTriangleCubaturePoints();
	/// constructor 
	HighOrderTriangleSetGeometryAlgorithms();
    virtual ~HighOrderTriangleSetGeometryAlgorithms() {}
public:
	virtual void init();
	virtual void reinit();
    virtual void draw(const core::visual::VisualParams* vparams);
	/// returns a pointer to the BezierTriangleSetTopologyContainer object
	virtual HighOrderTriangleSetTopologyContainer *getTopologyContainer() const {
		return container;
	}
	/// computes the nodal value given the tetrahedron index, the barycentric coordinates and the vector of nodal values
	virtual Coord computeNodalValue(const size_t triangleIndex,const Vec3 barycentricCoordinate,const VecCoord& p); 
	/// computes the nodal value assuming that the position is the regular position in the mechanical state object
	virtual Coord computeNodalValue(const size_t triangleIndex,const Vec3 barycentricCoordinate); 
	/// computes the shape function 
	virtual Real computeShapeFunction(const TriangleIndexVector tbi, const Vec3 barycentricCoordinate)=0;
	/// computes the shape function derivative
    virtual Vec3 computeShapeFunctionDerivatives(const TriangleIndexVector tbi, const Vec3 barycentricCoordinate)=0;
    /// computes the shape function hessian
    virtual Mat33 computeShapeFunctionHessian(const TriangleIndexVector tbi, const Vec3 barycentricCoordinate)=0;
	/// computes Jacobian i.e. determinant of dpos/dparam
	virtual Real computeJacobian(const size_t triangleIndex, const Vec3 barycentricCoordinate,const VecCoord& p);
	/// compute Jacobian i.e. determinant of dpos/dparam
	virtual Real computeJacobian(const size_t triangleIndex, const Vec3 barycentricCoordinate);
	/// compute the derivative of the nodal value along the 4 barycentric coordinates ie dpos/dparami 
	virtual void computeNodalValueDerivatives(const size_t triangleIndex, const Vec3 barycentricCoordinate, const VecCoord& p,Deriv point[3]);

	/// test if the Bezier tetrahedron is a simple affine tesselation of a regular tetrahedron
	bool isBezierTriangleAffine(const size_t triangleIndex,const VecCoord& p, const Real tolerance=(Real)1e-5) const; 
	/// returns the mass coefficient used to compute the mass matrix when a tetrahedron is affine 
	virtual Real getAffineMassCoefficient(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2)=0;
	/// returns the mass coefficient used to compute the exact mass matrix for  a regular high order element  
	virtual Real getRegularMassCoefficient(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2,const size_t indr[2])=0;
	/// returns the stiffness coefficient used to compute the stiffness matrix when a tetrahedron is affine 
	virtual Mat33 getAffineStiffnessCoefficientMatrix(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2)=0;
		
protected:
    Data<bool> drawControlPointsEdges;
	Data<bool> drawSmoothEdges;
	Data<bool> drawControlPoints;
	Data<Real> d_referenceRadius; // radius to draw control points
};


} // namespace topology

} // namespace component

} // namespace sofa

#endif
