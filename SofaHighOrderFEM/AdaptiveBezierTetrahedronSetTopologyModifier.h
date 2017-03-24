
#ifndef SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERTETRAHEDRONSETTOPOLOGYMODIFIER_H
#define SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERTETRAHEDRONSETTOPOLOGYMODIFIER_H


#include <SofaBaseTopology/TetrahedronSetTopologyModifier.h>
#include "initHighOrderFEM.h"
#include <sofa/helper/map.h>

namespace sofa
{

namespace component
{

namespace topology
{




/** indicates that some points are about to be removed */
class SOFA_HIGHORDER_FEM_API HighOrderTetrahedronDegreeChanged : public core::topology::TopologyChange
{
public:
    HighOrderTetrahedronDegreeChanged(const sofa::helper::vector<unsigned int>& _eArray,
		const sofa::helper::vector<unsigned int>& _tArray,
		const sofa::helper::vector<unsigned int>& _tetArray) : core::topology::TopologyChange(core::topology::TOPOLOGYCHANGE_LASTID),
        highOrderTetrahedronDegreeChangedEdgeArray(_eArray),
		highOrderTetrahedronDegreeChangedTrangleArray(_tArray),
		highOrderTetrahedronDegreeChangedTetrahedronArray(_tetArray)
    { }

	virtual ~HighOrderTetrahedronDegreeChanged(){}

	const sofa::helper::vector<unsigned int> &getEdgeArray() const { return highOrderTetrahedronDegreeChangedEdgeArray;	}
	const sofa::helper::vector<unsigned int> &getTriangleArray() const { return highOrderTetrahedronDegreeChangedTrangleArray;	}
    const sofa::helper::vector<unsigned int> &getTetrahedronArray() const { return highOrderTetrahedronDegreeChangedTetrahedronArray;	}


public:

	sofa::helper::vector<unsigned int> highOrderTetrahedronDegreeChangedEdgeArray;
	sofa::helper::vector<unsigned int>  highOrderTetrahedronDegreeChangedTrangleArray;
    sofa::helper::vector<unsigned int>  highOrderTetrahedronDegreeChangedTetrahedronArray;
};




/**
* A class that performs topology algorithms on an Adaptive Bezier Tetrahedron Set.
*/

class SOFA_HIGHORDER_FEM_API AdaptiveHighOrderTetrahedronSetTopologyModifier : public TetrahedronSetTopologyModifier
{
public:
	SOFA_CLASS(AdaptiveHighOrderTetrahedronSetTopologyModifier,TetrahedronSetTopologyModifier);
protected:
	AdaptiveHighOrderTetrahedronSetTopologyModifier()
		: TetrahedronSetTopologyModifier() 
	{}

	virtual ~AdaptiveHighOrderTetrahedronSetTopologyModifier() {}
public:
	virtual void init();

	virtual void reinit();

	/** \brief change the Bezier degree of an array of tetrahedra
	@param tetrahedra an array of vertex indices describing the tetrahedra to be created
	*/
	virtual void changeDegree(const sofa::helper::vector< unsigned int > &edges,
		const sofa::helper::vector< unsigned int > &triangles,
		const sofa::helper::vector< unsigned int > &tetrahedra);


};


} // namespace topology

} // namespace component

} // namespace sofa

#endif
