#ifndef SOFA_HIGHORDERTOPOLOGY_LAGRANGETETRAHEDRONSETGEOMETRYALGORITHMS_H
#define SOFA_HIGHORDERTOPOLOGY_LAGRANGETETRAHEDRONSETGEOMETRYALGORITHMS_H


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
* A class that provides geometry information on a high order Lagrange Tetrahedral Mesh
*/
template < class DataTypes >
class LagrangeTetrahedronSetGeometryAlgorithms : public HighOrderTetrahedronSetGeometryAlgorithms<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(LagrangeTetrahedronSetGeometryAlgorithms,DataTypes),SOFA_TEMPLATE(HighOrderTetrahedronSetGeometryAlgorithms,DataTypes));
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef sofa::defaulttype::Vec<4,Real> Vec4;
    typedef sofa::defaulttype::Mat<4,4,Real> Mat44;

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

	/// constructor 
	LagrangeTetrahedronSetGeometryAlgorithms();
    virtual ~LagrangeTetrahedronSetGeometryAlgorithms() {}
public:
	virtual void init();

	/// computes the shape function 
	virtual Real computeShapeFunction(const TetrahedronIndexVector tbi, const Vec4 barycentricCoordinate);
	/// computes the shape function gradient
    virtual Vec4 computeShapeFunctionDerivatives(const TetrahedronIndexVector tbi, const Vec4 barycentricCoordinate);
    /// computes the shape function hessian
    virtual Mat44 computeShapeFunctionHessian(const TetrahedronIndexVector tbi, const Vec4 barycentricCoordinate);

	/// returns the mass coefficient used to compute the mass matrix when a tetrahedron is affine 
	virtual Real getAffineMassCoefficient(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2);
	/// returns the mass coefficient used to compute the exact mass matrix for  a regular high order element  
	virtual Real getRegularMassCoefficient(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2,const size_t indr[3]);
	/// returns the stiffness coefficient used to compute the stiffness matrix when a tetrahedron is affine 
	virtual Mat44 getAffineStiffnessCoefficientMatrix(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2);
private:
	helper::vector< helper::vector<Real> > stirlingNumberArray;
	void checkCoefficients();	
	long int sterlingVector(const TetrahedronIndexVector tbiIn,const TetrahedronIndexVector tbiIk) const;
	/// compute  the mass coefficient used for  the exact mass matrix for  a regular high order element  
	 Real computeRegularMassCoefficient(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2,const size_t indr[3]);
	/// compute the mass coefficient for  the mass matrix when a tetrahedron is affine 
	 Real computeAffineMassCoefficient(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2);
	/// compute the stiffness coefficient for the stiffness matrix when a tetrahedron is affine 
	Mat44 computeAffineStiffnessCoefficientMatrix(const TetrahedronIndexVector tbi1, const TetrahedronIndexVector tbi2);

};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_HIGHORDERTOPOLOGY_LAGRANGETETRAHEDRONSETGEOMETRYALGORITHMS_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTetrahedronSetGeometryAlgorithms<defaulttype::Vec3dTypes>;
extern template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTetrahedronSetGeometryAlgorithms<defaulttype::Vec2dTypes>;
extern template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTetrahedronSetGeometryAlgorithms<defaulttype::Vec1dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTetrahedronSetGeometryAlgorithms<defaulttype::Vec3fTypes>;
extern template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTetrahedronSetGeometryAlgorithms<defaulttype::Vec2fTypes>;
extern template class SOFA_HIGHORDER_TOPOLOGY_API LagrangeTetrahedronSetGeometryAlgorithms<defaulttype::Vec1fTypes>;
#endif
#endif

} // namespace topology

} // namespace component

} // namespace sofa

#endif
