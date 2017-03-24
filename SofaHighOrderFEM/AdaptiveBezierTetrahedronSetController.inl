
#ifndef SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERTETRAHEDRONSETCONTTROLLER_H
#define SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERTETRAHEDRONSETCONTTROLLER_H


#include <sofa/core/behavior/BaseController.h>
#include "initHighOrderFEM.h"


namespace sofa
{

namespace component
{

namespace behavior
{
class AdaptiveBezierTetrahedronSetTopologyContainer;




/**
* A class that performs topology algorithms on an Adaptive Bezier Tetrahedron Set.
*/
template < class DataTypes >
class SOFA_HIGHORDER_FEM_API AdaptiveBezierTetrahedronSetController : public BaseController 
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(AdaptiveBezierTetrahedronSetController,DataTypes),
		BaseController);

	typedef typename DataTypes::Real Real;
	typedef typename DataTypes::Coord       Coord;
	typedef typename DataTypes::VecCoord VecCoord;
	typedef typename DataTypes::Deriv       Deriv;
    typedef typename sofa::defaulttype::Vec<4,Real> Vec4;
	typedef sofa::helper::vector<size_t> SeqBezierDegree;
	
protected:
    AdaptiveBezierTetrahedronSetController()
        : BaseController()
    {}

    virtual ~AdaptiveBezierTetrahedronSetController() {}
public:
    virtual void init();
	

protected:
	/** the object where the mechanical DOFs are stored */
	sofa::core::behavior::MechanicalState<DataTypes> *object;
    AdaptiveBezierTetrahedronSetTopologyContainer*					m_adaptiveContainer;
    BezierTetrahedronSetGeometryAlgorithms< DataTypes >*		m_bezierGeometryAlgorithms;

};


} // namespace topology

} // namespace component

} // namespace sofa

#endif
