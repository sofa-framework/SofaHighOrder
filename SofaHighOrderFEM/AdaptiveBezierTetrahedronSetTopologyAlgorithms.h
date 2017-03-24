
#ifndef SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERTETRAHEDRONSETTOPOLOGYALGORITHMS_H
#define SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERTETRAHEDRONSETTOPOLOGYALGORITHMS_H

#include "HighOrderTetrahedronSetTopologyContainer.h"
#include <SofaBaseTopology/TetrahedronSetTopologyAlgorithms.h>
#include "initHighOrderFEM.h"
#include <sofa/helper/map.h>

namespace sofa
{

namespace component
{

namespace topology
{
class AdaptiveHighOrderTetrahedronSetTopologyContainer;
class AdaptiveHighOrderTetrahedronSetTopologyModifier;



template < class DataTypes >
class BezierTetrahedronSetGeometryAlgorithms;



/**
* A class that performs topology algorithms on an Adaptive Bezier Tetrahedron Set.
*/
template < class DataTypes >
class SOFA_HIGHORDER_FEM_API AdaptiveBezierTetrahedronSetTopologyAlgorithms : public TetrahedronSetTopologyAlgorithms<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(AdaptiveBezierTetrahedronSetTopologyAlgorithms,DataTypes),
		SOFA_TEMPLATE(TetrahedronSetTopologyAlgorithms,DataTypes));

	typedef typename DataTypes::Real Real;
	typedef typename DataTypes::Coord       Coord;
	typedef typename DataTypes::VecCoord VecCoord;
	typedef typename DataTypes::Deriv       Deriv;
    
    typedef BaseMeshTopology::TetraID TetraID;
    typedef BaseMeshTopology::Tetra Tetra;
	typedef BaseMeshTopology::Edge Edge;
	typedef BaseMeshTopology::Triangle Triangle;
    typedef BaseMeshTopology::Tetrahedron Tetrahedron;
    typedef BaseMeshTopology::SeqTetrahedra SeqTetrahedra;
    typedef BaseMeshTopology::EdgesInTetrahedron EdgesInTetrahedron;
    typedef BaseMeshTopology::TrianglesInTetrahedron TrianglesInTetrahedron;
    typedef BaseMeshTopology::TetrahedraAroundTriangle TetrahedraAroundTriangle;
    typedef BaseMeshTopology::TetrahedraAroundEdge TetrahedraAroundEdge;

//    typedef typename DataTypes::Coord Coord;
	typedef sofa::helper::vector<size_t> SeqBezierDegree;
	typedef HighOrderTetrahedronSetTopologyContainer::HighOrderTetrahedronPointLocation HighOrderTetrahedronPointLocation;
	typedef std::pair<size_t,std::pair<HighOrderTetrahedronPointLocation,size_t> > ControlPointLocation;
	
protected:
    AdaptiveBezierTetrahedronSetTopologyAlgorithms();

    virtual ~AdaptiveBezierTetrahedronSetTopologyAlgorithms() {}
public:
    virtual void init();
	/** The tetrahedra whose ID are loweringDegreeTetrahedra will become linear while
	tetrahedra whose ID  raisingDegreeTetrahedra will become full degree */
	void updateTetrahedronDegree(sofa::helper::vector<TetraID>& loweringDegreeTetrahedra,
		sofa::helper::vector<TetraID>& raisingDegreeTetrahedra);


protected:
	/** the object where the mechanical DOFs are stored */
	sofa::core::behavior::MechanicalState<DataTypes> *object;
    AdaptiveHighOrderTetrahedronSetTopologyContainer*					m_adaptiveContainer;
	AdaptiveHighOrderTetrahedronSetTopologyModifier*					m_adaptiveModifier;
    BezierTetrahedronSetGeometryAlgorithms< DataTypes >*		m_bezierGeometryAlgorithms;
	std::vector<Coord > restPositionArray;
	std::map<ControlPointLocation,Coord > restPositionMap;
	std::map<ControlPointLocation,Coord  > offsetPositionMap;

};


} // namespace topology

} // namespace component

} // namespace sofa

#endif
