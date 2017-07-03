
#ifndef SOFA_HIGHORDERTOPOLOGY_HIGHORDERTETRAHEDRONSETTOPOLOGYCONTAINER_H
#define SOFA_HIGHORDERTOPOLOGY_HIGHORDERTETRAHEDRONSETTOPOLOGYCONTAINER_H
#include "initHighOrderTopology.h"

#include <SofaBaseTopology/TetrahedronSetTopologyContainer.h>
#include "HighOrderTopologyIndex.h"

namespace sofa
{

namespace component
{

namespace topology
{
class HighOrderTetrahedronSetTopologyModifier;

using core::topology::BaseMeshTopology;



/** a class that stores a set of high order  tetrahedra and provides access with adjacent triangles, edges and vertices 
A high order tetrahedron has exactly the same topology as a Tetrahedron but with additional (control) points on its edges, triangles and inside 
We use a Vec4D to number the control points inside  a high order tetrahedron */
class SOFA_HIGHORDER_TOPOLOGY_API HighOrderTetrahedronSetTopologyContainer : public TetrahedronSetTopologyContainer
{

public: 
    SOFA_CLASS(HighOrderTetrahedronSetTopologyContainer,TetrahedronSetTopologyContainer);

	typedef BaseMeshTopology::PointID			         PointID;
	typedef BaseMeshTopology::EdgeID			            EdgeID;
	typedef BaseMeshTopology::TriangleID		         TriangleID;
	typedef BaseMeshTopology::TetraID			         TetraID;
	typedef BaseMeshTopology::Edge				         Edge;
	typedef BaseMeshTopology::Triangle			         Triangle;
	typedef BaseMeshTopology::Tetra				         Tetra;
	typedef BaseMeshTopology::SeqTetrahedra			   SeqTetrahedra;




	typedef Tetra			Tetrahedron;
	typedef EdgesInTetrahedron		EdgesInTetrahedron;
	typedef TrianglesInTetrahedron	TrianglesInTetrahedron;
	typedef sofa::helper::vector<PointID>         VecPointID;
	typedef sofa::helper::vector<TetraID>         VecTetraID;
	typedef sofa::helper::vector<SReal> SeqWeights;
	typedef sofa::helper::vector<bool> SeqBools;
	typedef VecPointID  DOFInHighOrderTetrahedron;
	typedef sofa::helper::vector< VecPointID > SeqDOFInHighOrderTetrahedron;


	
	typedef sofa::defaulttype::Vec<4,int> ElementTetrahedronIndex;
	typedef sofa::defaulttype::Vec<4,size_t> LocalTetrahedronIndex;

	/* specify for each control point lying on an edge : the control point index, the index of the  edge, 
	the 2 integers specifying the position within this edge (i.e. 11 for a quadratic edge, 13 within a quartic edge).. */
	typedef sofa::helper::fixed_array<PointID,4> HighOrderEdgePosition;
	/* specify for each control point lying on a triangle  : the control point index, the index of the  triangle, 
	the 3 integers specifying the position within this triangle (i.e. 111 for a cubic triangle , 121 within a quartic triangle).. */
	typedef sofa::helper::fixed_array<PointID,5> HighOrderTrianglePosition;
	/* specify for each control point lying on a tetrahedron  : the control point index, the index of the  triangle, 
	the 4 integers specifying the position within this tetrahedron (i.e. 1111 for a quartic tetrahedron).. */
	typedef sofa::helper::fixed_array<PointID,6> HighOrderTetrahedronPosition;


	friend class HighOrderTetrahedronSetTopologyModifier;
	friend class Mesh2HighOrderTopologicalMapping;

public :
	// specifies where a high order control Point can lies with respect to the underlying tetrahedral mesh
	enum HighOrderTetrahedronPointLocation
    {
        POINT = 0,
        EDGE =1 ,
        TRIANGLE = 2,
        TETRAHEDRON = 3
    };
	
	/* describes where a control point is located : first integer is the index of the element on which the control point is lying,
	 then a pair indicating the type of element and its offset within this element */
	typedef std::pair<size_t,std::pair< HighOrderTetrahedronPointLocation,size_t> > ControlPointLocation;
protected:
    HighOrderTetrahedronSetTopologyContainer();

    virtual ~HighOrderTetrahedronSetTopologyContainer() {}
public:
    virtual void init();
	// build some maps specific of the degree of the tetrahedral elements.
	virtual void reinit();

    /// High order Specific Information Topology API
    /// @{
public :
	/// the degree of the high order Tetrahedron 1=linear, 2=quadratic...
	Data <size_t> d_degree;
	/// the number of control points corresponding to the vertices of the tetrahedra (different from the total number of points)
    Data<size_t> d_numberOfTetrahedralPoints;
	/// whether the high order tetrahedron is integral (false = classical high order splines) or rational splines (true)
	Data <SeqBools> d_isRationalSpline;
	/// the array of weights for rational splines
	Data <SeqWeights > d_weightArray;
	/// the list of edges used to describe the position of high order control points
    Data< helper::vector< Edge > > inputEdges;
	/// the list of triangles used to describe the position of high order control points
    Data< helper::vector< Triangle > > inputTriangles;
	/// the array describing the local position of a control point within an edge
    Data< helper::vector< HighOrderEdgePosition > > inputHighOrderEdgePositions;
	/// the array describing the local position of a control point within a triangle
    Data< helper::vector< HighOrderTrianglePosition > > inputHighOrderTrianglePositions;
	/// the array describing the local position of a control point within a tetrahedron
    Data< helper::vector< HighOrderTetrahedronPosition > > inputHighOrderTetrahedronPositions;
public :

	/// get the Degree of the high order Tetrahedron 
	HighOrderDegreeType getDegree() const;
	/// get the number of control points corresponding to the vertices of the tetrahedra 
	size_t getNumberOfTetrahedralPoints() const;
	/// get the global index of the control point associated with a given tetrahedron index and given its 4D Index   
	size_t getGlobalIndexOfControlPoint(const TetraID tetrahedronIndex,const TetrahedronIndexVector id) ;
    /// return   the global index of a control point from its location, the element index and its offset
    bool getGlobalIndexFromLocation(const HighOrderTetrahedronPointLocation location,
        const size_t elementIndex, const size_t elementOffset, size_t &globalIndex);
	/// get the indices of all control points associated with a given tetrahedron
	const VecPointID &getGlobalIndexArrayOfControlPoints(const TetraID tetrahedronIndex) const;
	/// return the tetrahedron index given the local index in a tetrahedron
	TetrahedronIndexVector getTetrahedronIndexVector(const size_t localIndex) const;
	/// get the Tetrahedron  Index Array of degree d
	sofa::helper::vector<TetrahedronIndexVector> getTetrahedronIndexArray() const;
	/// get the Tetrahedron  Index Array of a given degree 
	sofa::helper::vector<TetrahedronIndexVector> getTetrahedronIndexArrayOfGivenDegree(const HighOrderDegreeType deg) const;
	/** create an array which maps the local index of a Tetrahedron Index of degree d-1
	into a local index of a TBI of degree d by adding respectively (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1) **/
	sofa::helper::vector<LocalTetrahedronIndex> getMapOfTetrahedronIndexArrayFromInferiorDegree() const;
	/// return the local index in a tetrahedron from a tetrahedron index (inverse of getTetrahedronIndexVector())
	size_t getLocalIndexFromTetrahedronIndexVector(const TetrahedronIndexVector id) const;
	/// return the location, the element index and offset from the global index of a point
	void getLocationFromGlobalIndex(const size_t globalIndex, HighOrderTetrahedronPointLocation &location, 
		size_t &elementIndex, size_t &elementOffset) ;
		/// convert the edge offset into a EdgeIndex
	void getEdgeIndexVectorFromEdgeOffset(size_t offset, EdgeIndexVector &ebi);
	/// convert the triangle offset into a TriangleIndex
	void getTriangleIndexVectorFromTriangleOffset(size_t offset, TriangleIndexVector &tbi);
	/// convert the tetrahedron offset into a TriangleIndex
	void getTetrahedronIndexVectorFromTetrahedronOffset(size_t offset, TetrahedronIndexVector &tbi);
	/// check the control Point Topology
	virtual bool checkHighOrderTetrahedronTopology();
	/** \brief Returns the weight coordinate of the ith DOF. */
	virtual SReal getWeight(int i) const;
	/// returns the array of weights
	const SeqWeights & getWeightArray() const;
	// if the high order tetrahedron is rational or integral
	bool isRationalSpline(int i) const;
	// returns the Hierarchical index associated with the Lexicographic index i
	size_t getHierarchicalIndex(size_t i) const;
	// returns the Lexicographic index associated with the Hierarchical index i
	size_t getLexicographicIndex(size_t i) const;
	// returns the Lexicographic index associated with a Tetrahedron Vector index
	size_t getLexicographicIndex(const TetrahedronIndexVector tvi) const;
    /// Create element lists which are on topology border:
    virtual void createElementsOnBorder();
    /// @}


    inline friend std::ostream& operator<< (std::ostream& out, const HighOrderTetrahedronSetTopologyContainer& t)
    {
        helper::ReadAccessor< Data< sofa::helper::vector<Tetrahedron> > > m_tetrahedron = t.d_tetrahedron;
        out  << m_tetrahedron<< " "
                << t.m_edgesInTetrahedron<< " "
                << t.m_trianglesInTetrahedron;

        out << " "<< t.m_tetrahedraAroundVertex.size();
        for (unsigned int i=0; i<t.m_tetrahedraAroundVertex.size(); i++)
        {
            out << " " << t.m_tetrahedraAroundVertex[i];
        }
        out <<" "<< t.m_tetrahedraAroundEdge.size();
        for (unsigned int i=0; i<t.m_tetrahedraAroundEdge.size(); i++)
        {
            out << " " << t.m_tetrahedraAroundEdge[i];
        }
        out <<" "<< t.m_tetrahedraAroundTriangle.size();
        for (unsigned int i=0; i<t.m_tetrahedraAroundTriangle.size(); i++)
        {
            out << " " << t.m_tetrahedraAroundTriangle[i];
        }
		out << " " << t.d_degree.getValue();
		out << " " << t.getNbPoints();
        return out;
    }

    inline friend std::istream& operator>>(std::istream& in, HighOrderTetrahedronSetTopologyContainer& t)
    {
        unsigned int s=0;
        sofa::helper::vector< unsigned int > value;
        helper::WriteAccessor< Data< sofa::helper::vector<Tetrahedron> > > m_tetrahedron = t.d_tetrahedron;

        in >> m_tetrahedron >> t.m_edgesInTetrahedron >> t.m_trianglesInTetrahedron;


        in >> s;
        for (unsigned int i=0; i<s; i++)
        {
            in >> value;
            t.m_tetrahedraAroundVertex.push_back(value);
        }
        in >> s;
        for (unsigned int i=0; i<s; i++)
        {
            in >> value;
            t.m_tetrahedraAroundEdge.push_back(value);
        }
        in >> s;
        for (unsigned int i=0; i<s; i++)
        {
            in >> value;
            t.m_tetrahedraAroundTriangle.push_back(value);
        }
        HighOrderDegreeType bdt=0;
		in >> bdt;
		t.d_degree.setValue(bdt);
		int nbp;
		in >> nbp;
		t.setNbPoints(nbp);
        return in;
    }
protected:
	/// array describing the global  index of the DOFs used in weightedDOFArray
	SeqDOFInHighOrderTetrahedron  tetrahedronDOFArray;
	/// Map which provides the location (point, edge, triangle, tetrahedron) of a control point given its tetrahedron  index
	std::map<TetrahedronIndexVector,ElementTetrahedronIndex> elementMap;
	/// Map which provides the offset in the DOF vector for a control point lying on an edge 
	std::map<TetrahedronIndexVector,size_t> edgeOffsetMap;
	/// Map which provides the offset in the DOF vector for a control point lying on a triangle 
	std::map<TetrahedronIndexVector,size_t> triangleOffsetMap;
	/// Map which provides the offset in the DOF vector for a control point lying on a tetrahedron
	std::map<TetrahedronIndexVector,size_t> tetrahedronOffsetMap;
	/// Map which provides the rank in a control point from the array outputed by getGlobalIndexArrayOfControlPointsInTetrahedron (consistent with tetrahedronIndexArray) 
	std::map<TetrahedronIndexVector,size_t> localIndexMap;
	/// array of the tetrahedron  index outputed by the function getGlobalIndexArrayOfControlPointsInTetrahedron()
	sofa::helper::vector<TetrahedronIndexVector> tetrahedronIndexArray;
	/// array of the tetrahedron  index outputed by the function getGlobalIndexArrayOfControlPointsInTetrahedron()
	sofa::helper::vector<TetrahedronIndexVector> reducedDegreeTetrahedronIndexArray;
	/// convert triangle offset into triangle  index
	sofa::helper::vector<TriangleIndexVector> offsetToTriangleIndexArray;
	/// convert triangle offset into triangle  index
	sofa::helper::vector<TetrahedronIndexVector> offsetToTetrahedronIndexArray;
	/** Map which provides the global index of a control point knowing its location (i.e. triangle index and its TriangleIndex).
	This is empty by default since there is a default layout of control points based on edge and triangles indices */
	std::map<ControlPointLocation,size_t> locationToGlobalIndexMap;
	/** Map which provides the  location (i.e. triangle index and its TriangleIndexVector) of a control point knowing its  global index.
	Note that the location may not be unique.
	This is empty by default since there is a default layout of control points based on edge and triangles indices */
	std::map<size_t,ControlPointLocation> globalIndexToLocationMap;
	void getGlobalIndexArrayOfControlPointsInTetrahedron(const TetraID tetrahedronIndex, VecPointID & indexArray) ;
	// convert Hierarchical to Lexicographical order
	sofa::helper::vector<size_t> lexicographicIndexArray;
	// convert Lexicographical to Hierarchical  order
	sofa::helper::vector<size_t> hierarchicalIndexArray;
	// parse data from loader arrays
	void parseInputData();

};

} // namespace topology

} // namespace component

} // namespace sofa

#endif
