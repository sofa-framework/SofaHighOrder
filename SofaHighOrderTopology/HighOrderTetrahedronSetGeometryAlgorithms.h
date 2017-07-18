#ifndef SOFA_HIGHORDERTOPOLOGY_HIGHORDERTETRAHEDRONSETGEOMETRYALGORITHMS_H
#define SOFA_HIGHORDERTOPOLOGY_HIGHORDERTETRAHEDRONSETGEOMETRYALGORITHMS_H


#include <SofaBaseTopology/TetrahedronSetGeometryAlgorithms.h>
#include "HighOrderTetrahedronSetTopologyContainer.h"


namespace sofa
{

namespace component
{

namespace topology
{
using core::topology::BaseMeshTopology;


template< typename Real, int N> typename NumericalIntegrationDescriptor<Real, N>::QuadraturePointArray TetrahedronConicalRule(const size_t d);

/**
* A class that provides geometry information on an TetrahedronSet.
*/
template < class DataTypes >
class HighOrderTetrahedronSetGeometryAlgorithms : public TetrahedronSetGeometryAlgorithms<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(HighOrderTetrahedronSetGeometryAlgorithms,DataTypes),SOFA_TEMPLATE(TetrahedronSetGeometryAlgorithms,DataTypes));
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef sofa::defaulttype::Vec<4,Real> Vec4;
    typedef sofa::defaulttype::Mat<4,4,Real> Mat44;

	typedef BaseMeshTopology::Edge Edge;
	typedef BaseMeshTopology::PointID PointID;
	typedef BaseMeshTopology::TetraID TetraID;
	typedef BaseMeshTopology::Triangle Triangle;
	typedef BaseMeshTopology::Tetrahedron Tetrahedron;
	typedef BaseMeshTopology::SeqTetrahedra SeqTetrahedra;
	typedef BaseMeshTopology::SeqEdges SeqEdges;
	typedef BaseMeshTopology::TetrahedraAroundVertex TetrahedraAroundVertex;
	typedef BaseMeshTopology::TetrahedraAroundEdge TetrahedraAroundEdge;
	typedef BaseMeshTopology::TetrahedraAroundTriangle TetrahedraAroundTriangle;
	typedef BaseMeshTopology::EdgesInTetrahedron EdgesInTetrahedron;
	typedef BaseMeshTopology::TrianglesInTetrahedron TrianglesInTetrahedron;
	typedef HighOrderTetrahedronSetTopologyContainer::VecPointID VecPointID;

protected:
   
	/// container	
	HighOrderTetrahedronSetTopologyContainer *container; 
	/// degree of the polynomial
	HighOrderDegreeType degree; 
	// array of Tetrahedral Bezier indices
	sofa::helper::vector<TetrahedronIndexVector> tbiArray;
	/// the list of edges of the Bezier Tetrahedron used in the draw function
    std::set< std::pair<Edge,size_t> > bezierTetrahedronEdgeSet;
    bool initializedNewCubatureTables;
    void defineNewTetrahedronCubaturePoints();

	/// constructor 
	HighOrderTetrahedronSetGeometryAlgorithms();
    virtual ~HighOrderTetrahedronSetGeometryAlgorithms() {}
public:
	virtual void init();
	virtual void reinit();
    virtual void draw(const core::visual::VisualParams* vparams);
	/// returns a pointer to the BezierTetrahedronSetTopologyContainer object
	virtual HighOrderTetrahedronSetTopologyContainer *getTopologyContainer() const {
		return container;
	}
	/// computes the nodal value given the tetrahedron index, the barycentric coordinates and the vector of nodal values
	virtual Coord computeNodalValue(const size_t tetrahedronIndex,const Vec4 barycentricCoordinate,const VecCoord& p); 
	/// computes the nodal value assuming that the position is the regular position in the mechanical state object
	virtual Coord computeNodalValue(const size_t tetrahedronIndex,const Vec4 barycentricCoordinate); 
	/// computes the shape function 
	virtual Real computeShapeFunction(const TetrahedronIndexVector tbi, const Vec4 barycentricCoordinate)=0;
	/// computes the shape function derivative
    virtual Vec4 computeShapeFunctionDerivatives(const TetrahedronIndexVector tbi, const Vec4 barycentricCoordinate)=0;
    /// computes the shape function hessian
    virtual Mat44 computeShapeFunctionHessian(const TetrahedronIndexVector tbi, const Vec4 barycentricCoordinate)=0;
	/// computes Jacobian i.e. determinant of dpos/dparam
	virtual Real computeJacobian(const size_t tetrahedronIndex, const Vec4 barycentricCoordinate,const VecCoord& p);
	/// compute Jacobian i.e. determinant of dpos/dparam
	virtual Real computeJacobian(const size_t tetrahedronIndex, const Vec4 barycentricCoordinate);
	/// compute the derivative of the nodal value along the 4 barycentric coordinates ie dpos/dparami 
	virtual void computeNodalValueDerivatives(const size_t tetrahedronIndex, const Vec4 barycentricCoordinate, const VecCoord& p,Deriv point[4]);

	/// test if the Bezier tetrahedron is a simple affine tesselation of a regular tetrahedron
	bool isBezierTetrahedronAffine(const size_t tetrahedronIndex,const VecCoord& p, const Real tolerance=(Real)1e-5) const; 
	/// returns the mass coefficient used to compute the mass matrix when a tetrahedron is affine 
	virtual Real getAffineMassCoefficient(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2)=0;
	/// returns the mass coefficient used to compute the exact mass matrix for  a regular high order element  
	virtual Real getRegularMassCoefficient(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2,const size_t indr[3])=0;
	/// returns the stiffness coefficient used to compute the stiffness matrix when a tetrahedron is affine 
	virtual Mat44 getAffineStiffnessCoefficientMatrix(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2)=0;
		
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
