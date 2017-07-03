
#ifndef SOFA_HIGHORDERTOPOLOGY_BEZIERTRIANGLESETTOPOLOGYCONTAINER_H
#define SOFA_HIGHORDERTOPOLOGY_BEZIERTRIANGLESETTOPOLOGYCONTAINER_H

#include <SofaBaseTopology/TriangleSetTopologyContainer.h>
#include "HighOrderTopologyIndex.h"
#include "initHighOrderTopology.h"

namespace sofa
{

namespace component
{

namespace topology
{
class HighOrderTriangleSetTopologyModifier;

using core::topology::BaseMeshTopology;



typedef unsigned char HighOrderDegreeType;
typedef sofa::defaulttype::Vec<3,HighOrderDegreeType> TriangleIndexVector;


/** a class that stores a set of high order triangles and provides access with adjacent triangles, edges and vertices 
A high order triangle has exactly the same topology as a triangle  but with additional (control) points on its edges,  and inside 
We use a Vec3 to number the control points inside  a high order triangle  */
class SOFA_HIGHORDER_TOPOLOGY_API HighOrderTriangleSetTopologyContainer : public TriangleSetTopologyContainer
{
public:
	 SOFA_CLASS(HighOrderTriangleSetTopologyContainer,TriangleSetTopologyContainer);

	typedef BaseMeshTopology::PointID		            	PointID;
	typedef BaseMeshTopology::EdgeID		               	EdgeID;
	typedef BaseMeshTopology::TriangleID	               TriangleID;
	typedef BaseMeshTopology::Edge		        	         Edge;
	typedef BaseMeshTopology::Triangle	        	         Triangle;
	typedef BaseMeshTopology::SeqTriangles	        	      SeqTriangles;
	typedef BaseMeshTopology::EdgesInTriangle	         	EdgesInTriangle;
	typedef BaseMeshTopology::TrianglesAroundVertex    	    TrianglesAroundVertex;
	typedef BaseMeshTopology::TrianglesAroundEdge        	TrianglesAroundEdge;
	typedef sofa::helper::vector<TriangleID>                  VecTriangleID;
	typedef sofa::helper::vector<PointID>					  VecPointID;


	typedef sofa::defaulttype::Vec<3,int> ElementTriangleIndex;
	typedef sofa::defaulttype::Vec<3,size_t> LocalTriangleIndex;

	typedef sofa::helper::vector<SReal> SeqWeights;
	typedef sofa::helper::vector<bool> SeqBools;
	typedef sofa::helper::vector< VecPointID > SeqDOFInHighOrderTriangle;

	/* specify for each control point lying on an edge : the control point index, the index of the  edge, 
	 the 2 integers specifying the position within this edge (i.e. 11 for a quadratic edge, 13 within a quartic edge).. */
	typedef sofa::helper::fixed_array<PointID,4> HighOrderEdgePosition;
	/* specify for each control point lying on a triangle  : the control point index, the index of the  triangle, 
	 the 3 integers specifying the position within this triangle (i.e. 111 for a cubic triangle , 121 within a quartic triangle).. */
	typedef sofa::helper::fixed_array<PointID,5> HighOrderTrianglePosition;


    // specifies where a control Point can lies with respect to the underlying tetrahedral mesh
    enum HighOrderTrianglePointLocation
    {
        POINT = 0,
        EDGE = 1,
        TRIANGLE = 2,
        NONE = 3
    };
    //	typedef std::pair<size_t,TriangleIndexVector> ControlPointLocation;
    /* describes where a control point is located : first integer is the index of the element on which the control point is lying,
    then a pair indicating the type of element and its offset within this element */
    typedef std::pair<size_t, std::pair< HighOrderTrianglePointLocation, size_t> > ControlPointLocation;

    friend class HighOrderTriangleSetTopologyModifier;
	friend class Mesh2HighOrderTopologicalMapping;
	friend class HighOrderTetra2HighOrderTriangleTopologicalMapping;



protected:
    HighOrderTriangleSetTopologyContainer();

    virtual ~HighOrderTriangleSetTopologyContainer() {}
public:
    virtual void init();
	// build some maps specific of the degree of the tetrahedral elements.
	virtual void reinit();

    /// High order Specific Information Topology API
    /// @{
public :
	/// the degree of the high order triangle 1=linear, 2=quadratic...
	Data <size_t> d_degree;
	/// the number of control points corresponding to the vertices of the triangle mesh (different from the total number of points)
    Data<size_t> d_numberOfTriangularPoints;
	/// whether the high order triangles are integral (false = classical  splines) or rational splines (true)
	Data <SeqBools> d_isRationalSpline;
	/// the array of weights for rational splines
	Data <SeqWeights > d_weightArray;
	/// the list of edges used to describe the position of high order control points
    Data< helper::vector< Edge > > inputEdges;
	/// the array describing the local position of a control point within an edge
    Data< helper::vector< HighOrderEdgePosition > > inputHighOrderEdgePositions;
	/// the array describing the local position of a control point within a triangle
    Data< helper::vector< HighOrderTrianglePosition > > inputHighOrderTrianglePositions;
public :

	/// get the Degree of the high order triangle
	HighOrderDegreeType getDegree() const;
	/// get the number of control points corresponding to the vertices of the triangle mesh 
	size_t getNumberOfTriangularPoints() const;
	/// get the global index of the control point associated with a given tetrahedron index and given its 4D Index   
	size_t getGlobalIndexOfControlPoint(const TriangleID triangleIndex,const TriangleIndexVector id) ;
    /// return   the global index of a control point from its location, the element index and its offset
    bool getGlobalIndexFromLocation(const HighOrderTrianglePointLocation location,
        const size_t elementIndex, const size_t elementOffset, size_t &globalIndex);
	/// get the indices of all control points associated with a given triangle
	const VecPointID &getGlobalIndexArrayOfControlPoints(const TetraID triangleIndex) const;
	/// return the triangle index vector given the local index in a triangle
	TriangleIndexVector getTriangleIndexVector(const size_t localIndex) const;
	/// get the Triangle Index Array of degree d
	sofa::helper::vector<TriangleIndexVector> getTriangleIndexArray() const;
	/// get the Triangle  Index Array of a given degree 
	sofa::helper::vector<TriangleIndexVector> getTriangleIndexArrayOfGivenDegree(const HighOrderDegreeType deg) const;
	/** create an array which maps the local index of a Triangle Index of degree d-1
	into a local index of a TBI of degree d by adding respectively (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1) **/
	sofa::helper::vector<LocalTriangleIndex> getMapOfTriangleIndexArrayFromInferiorDegree() const;
	/** return the array describing each of the  (degree+1)*(degree+1) subtriangles with local indices ( i.e. indices between 0 and (degree+1)*(degree+2)/2  ) */ 
	sofa::helper::vector<LocalTriangleIndex> getLocalIndexSubtriangleArray() const;
	/** return the array describing each of the  (deg+1)*(deg+1) subtriangles with local indices ( i.e. indices between 0 and (deg+1)*(deg+2)/2  ) */ 
	sofa::helper::vector<LocalTriangleIndex> getLocalIndexSubtriangleArrayOfGivenDegree(const HighOrderDegreeType deg)  const;
	/// return the local index in a triangle from a triangle index (inverse of getTriangleIndex())
	size_t getLocalIndexFromTriangleIndexVector(const TriangleIndexVector id) const;
	/// return the location, the element index and offset from the global index of a point
	void getLocationFromGlobalIndex(const size_t globalIndex, HighOrderTrianglePointLocation &location, 
		size_t &elementIndex, size_t &elementOffset) ;
	/// convert the edge offset into a EdgeIndexVector
	void getEdgeIndexVectorFromEdgeOffset(size_t offset, EdgeIndexVector &ebi);
	/// convert the triangle offset into a TriangleIndexVector
	void getTriangleIndexVectorFromTriangleOffset(size_t offset, TriangleIndexVector &tbi);
	/// check the high order triangle Topology
	bool checkHighOrderTriangleTopology();
    /** \brief Returns the weight coordinate of the ith DOF. */
    virtual SReal getWeight(int i) const;
    /// returns the array of weights
    const SeqWeights & getWeightArray() const;
    // if the triangle is rational or integral
    bool isRationalSpline(int i) const;
    // returns the Hierarchical index associated with the Lexicographic index i
    size_t getHierarchicalIndex(size_t i) const;
    // returns the Lexicographic index associated with the Hierarchical index i
    size_t getLexicographicIndex(size_t i) const;
    // returns the Lexicographic index associated with a Triangle Vector index
    size_t getLexicographicIndex(const TriangleIndexVector tvi) const;
    /// Create element lists which are on topology border:
    virtual void createElementsOnBorder();
	 /// @}


protected:
	/// array describing the global  index of the DOFs used in weightedDOFArray
	SeqDOFInHighOrderTriangle  triangleDOFArray;
	/** Map which provides the global index of a control point knowing its location (i.e. triangle index and its TriangleIndexVector).
	This is empty by default since there is a default layout of control points based on edge and triangles indices */
	std::map<ControlPointLocation,size_t> locationToGlobalIndexMap;
	/** Map which provides the  location (i.e. triangle index and its TriangleIndexVector) of a control point knowing its  global index.
	Note that the location may not be unique.
	This is empty by default since there is a default layout of control points based on edge and triangles indices */
    std::map<size_t, ControlPointLocation> globalIndexToLocationMap;


	/// Map which provides the location (point, edge, Triangle) of a control point given its Triangle Bezier index
	std::map<TriangleIndexVector,ElementTriangleIndex> elementMap;
	/// Map which provides the offset in the DOF vector for a control point lying on an edge 
	std::map<TriangleIndexVector,size_t> edgeOffsetMap;
	/// Map which provides the offset in the DOF vector for a control point lying on a triangle 
	std::map<TriangleIndexVector,size_t> triangleOffsetMap;

	/// Map which provides the rank in a control point from the array outputted by getGlobalIndexArrayOfControlPointsInTriangle (consistent with bezierIndexArray) 
	std::map<TriangleIndexVector,size_t> localIndexMap;
	/// array of the Triangle index outputed by the function getGlobalIndexArrayOfControlPointsInTriangle()
	sofa::helper::vector<TriangleIndexVector> bezierIndexArray;
	/// array of the Triangle index outputed by the function getGlobalIndexArrayOfControlPointsInTriangle()
	sofa::helper::vector<TriangleIndexVector> reducedDegreeIndexArray;
	/// convert triangle offset into triangle bezier index
	sofa::helper::vector<TriangleIndexVector> offsetToTriangleIndexVectorArray;
	// convert Hierarchical to Lexicographical order
	sofa::helper::vector<size_t> lexicographicIndexArray;
	// convert Lexicographical to Hierarchical  order
	sofa::helper::vector<size_t> hierarchicalIndexArray;
	void getGlobalIndexArrayOfControlPointsInTriangle(const TriangleID triangleIndex, VecPointID & indexArray) ;
	void parseInputData(); 

};

} // namespace topology

} // namespace component

} // namespace sofa

#endif
