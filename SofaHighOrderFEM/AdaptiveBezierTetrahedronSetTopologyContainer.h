#ifndef SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERTETRAHEDRONSETTOPOLOGYCONTAINER_H
#define SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERTETRAHEDRONSETTOPOLOGYCONTAINER_H
//#include "initHighOrderFEM.h"

#include "HighOrderTetrahedronSetTopologyContainer.h"
#include "initHighOrderFEM.h"

namespace sofa
{

namespace component
{

namespace topology
{
class AdaptiveHighOrderTetrahedronSetTopologyModifier;
class AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping;

template < class DataTypes >
class AdaptiveBezierTetrahedronSetTopologyAlgorithms;

using core::topology::BaseMeshTopology;


/** a class that stores a set of Bezier tetrahedra and provides access with adjacent triangles, edges and vertices 
A Bezier Tetrahedron has exactly the same topology as a Tetrahedron but with additional (control) points on its edges, triangles and inside 
We use a Vec4D to number the control points inside  a Bezier tetrahedron */
class SOFA_HIGHORDER_FEM_API AdaptiveHighOrderTetrahedronSetTopologyContainer : public HighOrderTetrahedronSetTopologyContainer
{

public: 
    SOFA_CLASS(AdaptiveHighOrderTetrahedronSetTopologyContainer,HighOrderTetrahedronSetTopologyContainer);

	typedef BaseMeshTopology::PointID			         PointID;
	typedef BaseMeshTopology::EdgeID			            EdgeID;
	typedef BaseMeshTopology::TriangleID		         TriangleID;
	typedef BaseMeshTopology::TetraID			         TetraID;
	typedef BaseMeshTopology::Edge				         Edge;
	typedef BaseMeshTopology::Triangle			         Triangle;
	typedef BaseMeshTopology::Tetra				         Tetra;
	typedef BaseMeshTopology::SeqTetrahedra			   SeqTetrahedra;


	typedef sofa::defaulttype::Vec<4,SReal> Vec4;
	typedef Tetra			Tetrahedron;
	typedef EdgesInTetrahedron		EdgesInTetrahedron;
	typedef TrianglesInTetrahedron	TrianglesInTetrahedron;
	typedef sofa::helper::vector<PointID>         VecPointID;
	typedef sofa::helper::vector<TetraID>         VecTetraID;
	typedef sofa::helper::vector<SReal> SeqWeights;
	typedef std::pair<size_t,SReal>  WeightedDOF;
	typedef sofa::helper::vector<std::pair<size_t,SReal> > WeightedDOFArray;
	typedef sofa::helper::vector<bool> SeqBools;
	typedef sofa::helper::vector<size_t> SeqBezierDegree;
	typedef VecPointID  BezierDOFInTetrahedron;
	typedef sofa::helper::vector< VecPointID > SeqBezierDOFInTetrahedron;
	typedef sofa::helper::vector<WeightedDOFArray> SeqWeightedDOFArray;
	typedef sofa::helper::vector<sofa::helper::vector<WeightedDOFArray> > SeqSeqWeightedDOFArray;
	
	typedef HighOrderTetrahedronSetTopologyContainer::HighOrderTetrahedronPointLocation BezierTetrahedronPointLocation;
	typedef HighOrderTetrahedronSetTopologyContainer::ControlPointLocation  ControlPointLocation;

	friend class AdaptiveHighOrderTetrahedronSetTopologyModifier;
	friend class AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping;

	template < class DataTypes >
	friend class AdaptiveBezierTetrahedronSetTopologyAlgorithms;

protected:
    AdaptiveHighOrderTetrahedronSetTopologyContainer();

    virtual ~AdaptiveHighOrderTetrahedronSetTopologyContainer() {}
public:
    virtual void init();
	// build some maps specific of the degree of the tetrahedral elements.
	virtual void reinit();
	virtual bool checkHighOrderTetrahedronTopology();
	

    /// Adaptive Bezier Specific Information Topology API
    /// @{
public :
	/// the degree of the Bezier tetrahedra either 1=linear or degree
	Data <SeqBezierDegree> d_tetrahedronDegree;
		/// the degree of the Bezier triangles either 1=linear or degree
	Data <SeqBezierDegree> d_triangleDegree;
		/// the degree of the Bezier edges  either 1=linear or degree
	Data <SeqBezierDegree> d_edgeDegree;
	/// indicates is surface or volumetric extrapolation are used to update surface triangles
	Data<bool> d_useSurfaceExtrapolation;
public :
	/** In swda get  the local indices of the degrees of freedom and their weight that are mapped for each control 
	points of a Bezier tetrahedron. In da get the global DOF id of each local DOF index */ 
	void getWeightedDOFOfTetrahedron(const size_t tetrahedronIndex, SeqWeightedDOFArray &swda,BezierDOFInTetrahedron &);
	// returns a pointer to an array describing the DOFs of all tetrahedra
	const SeqSeqWeightedDOFArray &getWeightedDOFArrayArray() const;
	// returns a pointer to an array describing the DOFs of all tetrahedra
	const SeqBezierDOFInTetrahedron &getTetrahedronDOFArrayArray() const;
	// returns a pointer to an array describing the DOFs of a given tetrahedron
	const SeqWeightedDOFArray &getWeightedDOFArray(const unsigned int i);
	// returns a pointer to an array describing the DOFs of a given tetrahedron
	const BezierDOFInTetrahedron &getTetrahedronDOFArray(const unsigned int i);
	/// get the Degree of the Bezier Tetrahedron 
	size_t getTetrahedronDegree(const size_t tetrahedronIndex) const;
	/// get the array descriving the Degree of all Bezier Tetrahedra 
	const  SeqBezierDegree &getTetrahedronDegreeArray() const;
	/// get the Degree of the Bezier Tetrahedron 
	size_t getTriangleDegree(const size_t triangleIndex) const;
	/// get the array descriving the Degree of all Bezier triangles 
	const  SeqBezierDegree &getTriangleDegreeArray() const;
	/// get the Degree of the Bezier edge
	size_t getEdgeDegree(const size_t edgeIndex) const;
	/// get the array descriving the Degree of all Bezier edges 
	const  SeqBezierDegree &getEdgeDegreeArray() const;

	 /// @}
	// sets of triangles on surface
	std::set<size_t> triangleOnSurfaceSet;
	// set of edges on surface
	std::set<size_t> edgeOnSurfaceSet;
	std::map<size_t,size_t> trianglesOnSurfaceMap;
	std::map<size_t, std::pair<size_t,size_t> > edgeOnSurfaceMap;
	std::map<std::pair<size_t,size_t>,Vec4 > triangleControlPointsBarycentricCoord;
	std::map<std::pair<size_t,size_t>,std::pair<Vec4,Vec4> > edgeControlPointsBarycentricCoord;
   
protected:

	/// array describing the DOF involved in each tetrahedron and their weights : use local indices
	SeqSeqWeightedDOFArray  weightedDOFArray;
	// store the bezier weights to restore their value when a tetrahedron is degree raised
	std::map<ControlPointLocation,SReal> bezierWeightsMap;

};

} // namespace topology

} // namespace component

} // namespace sofa

#endif
