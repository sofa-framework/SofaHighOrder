
#ifndef SOFA_HIGHORDERTOPOLOGY_BEZIERTETRAHEDRONSETGEOMETRYALGORITHMS_H
#define SOFA_HIGHORDERTOPOLOGY_BEZIERTETRAHEDRONSETGEOMETRYALGORITHMS_H


#include <SofaBaseTopology/TetrahedronSetGeometryAlgorithms.h>
#include "HighOrderTetrahedronSetTopologyContainer.h"
#include "HighOrderTetrahedronSetGeometryAlgorithms.h"

namespace sofa
{

namespace component
{

namespace topology
{
using core::topology::BaseMeshTopology;


/**
* A class that provides geometry information on an TetrahedronSet.
*/
template < class DataTypes >
class BezierTetrahedronSetGeometryAlgorithms : public HighOrderTetrahedronSetGeometryAlgorithms<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(BezierTetrahedronSetGeometryAlgorithms,DataTypes),SOFA_TEMPLATE(HighOrderTetrahedronSetGeometryAlgorithms,DataTypes));
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
	typedef typename DataTypes::Deriv Deriv;
    typedef sofa::defaulttype::Vec<4,Real> Vec4;
    typedef sofa::defaulttype::Mat<4,4,Real> Mat44;

	typedef BaseMeshTopology::Edge Edge;
	typedef BaseMeshTopology::PointID PointID;
	typedef BaseMeshTopology::TetraID TetraID;
	typedef BaseMeshTopology::Tetra Tetra;
	typedef BaseMeshTopology::SeqTetrahedra SeqTetrahedra;
	typedef BaseMeshTopology::SeqEdges SeqEdges;
	typedef BaseMeshTopology::TetrahedraAroundVertex TetrahedraAroundVertex;
	typedef BaseMeshTopology::TetrahedraAroundEdge TetrahedraAroundEdge;
	typedef BaseMeshTopology::TetrahedraAroundTriangle TetrahedraAroundTriangle;
	typedef BaseMeshTopology::EdgesInTetrahedron EdgesInTetrahedron;
	typedef BaseMeshTopology::TrianglesInTetrahedron TrianglesInTetrahedron;
	typedef HighOrderTetrahedronSetTopologyContainer::VecPointID VecPointID;

	typedef Tetra Tetrahedron;
protected:
   
	// array of Bernstein coefficient following the same order as tbiArray
	sofa::helper::vector<Real> bernsteinCoefficientArray;
	// map used to store the Bernstein coefficient given a Tetrahedron Bezier Index
	std::map<TetrahedronIndexVector,Real> bernsteinCoeffMap;

	/// constructor 
	BezierTetrahedronSetGeometryAlgorithms();
    virtual ~BezierTetrahedronSetGeometryAlgorithms() {}
public:
	virtual void init();

	/// computes the shape function 
	virtual Real computeShapeFunction(const TetrahedronIndexVector tbi, const Vec4 barycentricCoordinate);
	/// computes the shape function gradient
    virtual Vec4 computeShapeFunctionDerivatives(const TetrahedronIndexVector tbi, const Vec4 barycentricCoordinate);
    /// computes the shape function hessian
    virtual Mat44 computeShapeFunctionHessian(const TetrahedronIndexVector tbi, const Vec4 barycentricCoordinate);
	/// computes Jacobian i.e. determinant of dpos/dmu
	virtual Real computeJacobian(const size_t tetrahedronIndex, const Vec4 barycentricCoordinate,const VecCoord& p);
	/// compute the derivative of the nodal value along the 4 barycentric coordinates ie dpos/dparami 
	virtual void computeNodalValueDerivatives(const size_t tetrahedronIndex, const Vec4 barycentricCoordinate, const VecCoord& p,Deriv point[4]);
	/// compute the 4 De Casteljeau  of degree d-1
	virtual void computeDeCasteljeauPoints(const size_t tetrahedronIndex, const Vec4 barycentricCoordinate, const VecCoord& p,Coord point[4]);
	/// returns the mass coefficient used to compute the mass matrix when a tetrahedron is affine 
	virtual Real getAffineMassCoefficient(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2);
	/// returns the mass coefficient used to compute the exact mass matrix for  a regular high order element  
	virtual Real getRegularMassCoefficient(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2,const size_t indr[3]);
	/// returns the stiffness coefficient used to compute the stiffness matrix when a tetrahedron is affine 
	virtual Mat44 getAffineStiffnessCoefficientMatrix(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2);
	   
    /// computes the shape function for any degree
    Real computeShapeFunctionOfGivenDegree(const TetrahedronIndexVector tbi, const Vec4 barycentricCoordinate, const HighOrderDegreeType deg);
    /// computes the shape function gradient for any degree
    Vec4 computeShapeFunctionDerivativesOfGivenDegree(const TetrahedronIndexVector tbi, const Vec4 barycentricCoordinate, const HighOrderDegreeType deg);
    template<class T>
    static bool canCreate(T*& obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
        return core::objectmodel::BaseObject::canCreate(obj, context, arg);
    }
private:

	void checkCoefficients();	

	/// compute  the mass coefficient used for  the exact mass matrix for  a regular high order element  
	Real computeRegularMassCoefficient(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2,const size_t indr[3]);
	/// compute the mass coefficient for  the mass matrix when a tetrahedron is affine 
	Real computeAffineMassCoefficient(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2);
	/// compute the stiffness coefficient for the stiffness matrix when a tetrahedron is affine 
	Mat44 computeAffineStiffnessCoefficientMatrix(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2);
	
};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_HIGHORDERTOPOLOGY_BEZIERTETRAHEDRONSETGEOMETRYALGORITHMS_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_HIGHORDER_TOPOLOGY_API BezierTetrahedronSetGeometryAlgorithms<defaulttype::Vec3dTypes>;
extern template class SOFA_HIGHORDER_TOPOLOGY_API BezierTetrahedronSetGeometryAlgorithms<defaulttype::Vec2dTypes>;
extern template class SOFA_HIGHORDER_TOPOLOGY_API BezierTetrahedronSetGeometryAlgorithms<defaulttype::Vec1dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_HIGHORDER_TOPOLOGY_API BezierTetrahedronSetGeometryAlgorithms<defaulttype::Vec3fTypes>;
extern template class SOFA_HIGHORDER_TOPOLOGY_API BezierTetrahedronSetGeometryAlgorithms<defaulttype::Vec2fTypes>;
extern template class SOFA_HIGHORDER_TOPOLOGY_API BezierTetrahedronSetGeometryAlgorithms<defaulttype::Vec1fTypes>;
#endif
#endif

} // namespace topology

} // namespace component

} // namespace sofa

#endif
