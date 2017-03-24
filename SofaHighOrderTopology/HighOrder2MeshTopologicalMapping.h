#ifndef SOFA_HIGHORDERTOPOLOGY_BEZIER2MESHTOPOLOGICALMAPPING_H
#define SOFA_HIGHORDERTOPOLOGY_BEZIER2MESHTOPOLOGICALMAPPING_H
#include "initHighOrderTopology.h"


#include <sofa/core/topology/TopologicalMapping.h>
#include <sofa/core/topology/Topology.h>

#include <sofa/defaulttype/Vec.h>
#include <map>
#include <set>



namespace sofa { namespace component { namespace mapping { template<typename  D, typename E> class Bezier2MeshMechanicalMapping; } } }


namespace sofa
{
namespace component
{
namespace topology
{
/**
 * This class, called HighOrder2MeshTopologicalMapping, is a specific implementation of the interface TopologicalMapping where :
 *
 * INPUT TOPOLOGY = A BezierTetrahedronSetTopology or BezierTriangleSetTopology as a tesselated version of the input mesh 
  * OUTPUT TOPOLOGY = a Tetrahedral or triangular mesh interpolated with a given degree of tesselation from its Bezier mesh
 *
 * This Topological mapping handles the specic input topology of Bezier elements and is made more efficient by using precomputations of maps
 *
 * HighOrder2MeshTopologicalMapping class is templated by the pair (INPUT TOPOLOGY, OUTPUT TOPOLOGY)
 *
*/

    class SOFA_HIGHORDER_TOPOLOGY_API HighOrder2MeshTopologicalMapping : public sofa::core::topology::TopologicalMapping
{
public:
    SOFA_CLASS(HighOrder2MeshTopologicalMapping,sofa::core::topology::TopologicalMapping);
	template<typename D, typename E>  friend class sofa::component::mapping::Bezier2MeshMechanicalMapping;
	typedef sofa::core::topology::Topology::Tetrahedron Tetrahedron;
	typedef sofa::core::topology::Topology::TetraID TetraID;
	typedef sofa::core::topology::Topology::Tetra Tetra;
	typedef sofa::core::topology::Topology::Edge Edge;
	typedef sofa::core::topology::Topology::Triangle Triangle;
	typedef sofa::core::topology::BaseMeshTopology::SeqTriangles SeqTriangles;
	typedef sofa::core::topology::BaseMeshTopology::EdgesInTriangle EdgesInTriangle;
	typedef sofa::core::topology::BaseMeshTopology::TrianglesAroundEdge TrianglesAroundEdge;

protected:
    /** \brief Constructor.
     *
     */
    HighOrder2MeshTopologicalMapping ();

    /** \brief Destructor.
     *
         * Does nothing.
         */
    virtual ~HighOrder2MeshTopologicalMapping();
public:
    /** \brief Initializes the target BaseTopology from the source BaseTopology.
     */
    virtual void init();

    /// create a number of subtetrahedra depending on the level of tesselation. Degree is the number of times an edge will be split. 
	Data < unsigned int > d_tesselationTetrahedronDegree;
	/// create a number of subtriangles depending on the level of tesselation. Degree is the number of times an edge will be split. 
	Data < unsigned int > d_tesselationTriangleDegree;
protected:
	// local indexing of points inside tessellated triangles
	sofa::helper::vector<sofa::defaulttype::Vec<3,unsigned char > > tesselatedTriangleIndices; 
	sofa::helper::vector<sofa::core::topology::Topology::Edge > edgeTriangleArray; 
	sofa::helper::vector<sofa::helper::vector<size_t> > bezierEdgeArray; 
	/// for each macro triangle set the index of tesselated points inside that triangle (used for nmal computation)
	sofa::helper::vector< sofa::helper::vector<size_t> > globalIndexTesselatedBezierTriangleArray; 
	sofa::helper::vector<size_t> local2GlobalBezierVertexArray;
	sofa::helper::vector<int> global2LocalBezierVertexArray;
	// the number of points in the output triangulation
	size_t nbPoints;
public :
	/// Method called at each topological changes propagation which comes from the INPUT topology to adapt the OUTPUT topology :
	virtual void updateTopologicalMappingTopDown();

};



} // namespace topology
} // namespace component
} // namespace sofa

#endif // SOFA_HIGHORDERTOPOLOGY_MESH2BEZIERTOPOLOGICALMAPPING_H
