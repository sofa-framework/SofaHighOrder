#ifndef SOFA_HIGHORDERTOPOLOGY_BEZIERTRIANGLESETGEOMETRYALGORITHMS_H
#define SOFA_HIGHORDERTOPOLOGY_BEZIERTRIANGLESETGEOMETRYALGORITHMS_H


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
* A class that provides geometry information on an TriangleSet.
*/
template < class DataTypes >
class BezierTriangleSetGeometryAlgorithms : public HighOrderTriangleSetGeometryAlgorithms<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(BezierTriangleSetGeometryAlgorithms,DataTypes),SOFA_TEMPLATE(HighOrderTriangleSetGeometryAlgorithms,DataTypes));
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
	typedef BaseMeshTopology::SeqTriangles SeqTriangles;
	typedef BaseMeshTopology::SeqEdges SeqEdges;
	typedef BaseMeshTopology::TrianglesAroundVertex TrianglesAroundVertex;
	typedef BaseMeshTopology::TrianglesAroundEdge TrianglesAroundEdge;
	typedef BaseMeshTopology::EdgesInTriangle EdgesInTriangle;
	typedef HighOrderTriangleSetTopologyContainer::VecPointID VecPointID;

protected:
   
	// array of Bernstein coefficient following the same order as tbiArray
	sofa::helper::vector<Real> bernsteinCoefficientArray;
	// map used to store the Bernstein coefficient given a Triangle Bezier Index
	std::map<TriangleIndexVector,Real> bernsteinCoeffMap;

	/// constructor 
	BezierTriangleSetGeometryAlgorithms();
    virtual ~BezierTriangleSetGeometryAlgorithms() {}
public:
	virtual void init();

	/// computes the shape function 
	virtual Real computeShapeFunction(const TriangleIndexVector tbi, const Vec3 barycentricCoordinate);
	/// computes the shape function gradient
    virtual Vec3 computeShapeFunctionDerivatives(const TriangleIndexVector tbi, const Vec3 barycentricCoordinate);
    /// computes the shape function hessian
    virtual Mat33 computeShapeFunctionHessian(const TriangleIndexVector tbi, const Vec3 barycentricCoordinate);
	/// computes Jacobian i.e. determinant of dpos/dmu
	virtual Real computeJacobian(const size_t triangleIndex, const Vec3 barycentricCoordinate,const VecCoord& p);
	/// compute the derivative of the nodal value along the 4 barycentric coordinates ie dpos/dparami 
	virtual void computeNodalValueDerivatives(const size_t triangleIndex, const Vec3 barycentricCoordinate, const VecCoord& p,Deriv point[4]);
	/// compute the 4 De Casteljeau  of degree d-1
	virtual void computeDeCasteljeauPoints(const size_t triangleIndex, const Vec3 barycentricCoordinate, const VecCoord& p,Coord point[4]);
	/// returns the mass coefficient used to compute the mass matrix when a triangle is affine 
	virtual Real getAffineMassCoefficient(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2);
	/// returns the mass coefficient used to compute the exact mass matrix for  a regular high order element  
	virtual Real getRegularMassCoefficient(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2,const size_t indr[2]);
	/// returns the stiffness coefficient used to compute the stiffness matrix when a triangle is affine 
	virtual Mat33 getAffineStiffnessCoefficientMatrix(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2);
	    template<class T>
    static bool canCreate(T*& obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
        return core::objectmodel::BaseObject::canCreate(obj, context, arg);
    }
private:

	void checkCoefficients();	

	/// compute  the mass coefficient used for  the exact mass matrix for  a regular high order element  
	Real computeRegularMassCoefficient(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2,const size_t indr[2]);
	/// compute the mass coefficient for  the mass matrix when a triangle is affine 
	Real computeAffineMassCoefficient(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2);
	/// compute the stiffness coefficient for the stiffness matrix when a triangle is affine 
	Mat33 computeAffineStiffnessCoefficientMatrix(const TriangleIndexVector tbi1, const TriangleIndexVector tbi2);
	
};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_HIGHORDERTOPOLOGY_BEZIERTRIANGLESETGEOMETRYALGORITHMS_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_HIGHORDER_TOPOLOGY_API BezierTriangleSetGeometryAlgorithms<defaulttype::Vec3dTypes>;
extern template class SOFA_HIGHORDER_TOPOLOGY_API BezierTriangleSetGeometryAlgorithms<defaulttype::Vec2dTypes>;
extern template class SOFA_HIGHORDER_TOPOLOGY_API BezierTriangleSetGeometryAlgorithms<defaulttype::Vec1dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_HIGHORDER_TOPOLOGY_API BezierTriangleSetGeometryAlgorithms<defaulttype::Vec3fTypes>;
extern template class SOFA_HIGHORDER_TOPOLOGY_API BezierTriangleSetGeometryAlgorithms<defaulttype::Vec2fTypes>;
extern template class SOFA_HIGHORDER_TOPOLOGY_API BezierTriangleSetGeometryAlgorithms<defaulttype::Vec1fTypes>;
#endif
#endif

} // namespace topology

} // namespace component

} // namespace sofa

#endif
