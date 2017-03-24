
#ifndef SOFA_COMPONENT_TOPOLOGY_ADAPTIVEBEZIERTETRA2BEZIERTRIANGLETOPOLOGICALMAPPING_H
#define SOFA_COMPONENT_TOPOLOGY_ADAPTIVEBEZIERTETRA2BEZIERTRIANGLETOPOLOGICALMAPPING_H

#include <sofa/core/topology/TopologicalMapping.h>
#include "initHighOrderFEM.h"
#include <sofa/core/topology/BaseMeshTopology.h>

namespace sofa
{


namespace component
{
	namespace mapping { template<typename  D, typename E> class AdaptiveBezierTetra2BezierTriangleMechanicalMapping; } 
namespace topology
{

	class AdaptiveHighOrderTetrahedronSetTopologyContainer;
	class HighOrderTriangleSetTopologyContainer;
/**
 * This class, called Tetra2TriangleTopologicalMapping, is a specific implementation of the interface TopologicalMapping where :
 *
 * INPUT TOPOLOGY = TetrahedronSetTopology
 * OUTPUT TOPOLOGY = TriangleSetTopology, as the boundary of the INPUT TOPOLOGY
 *
 * Tetra2TriangleTopologicalMapping class is templated by the pair (INPUT TOPOLOGY, OUTPUT TOPOLOGY)
 *
*/

class SOFA_HIGHORDER_FEM_API AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping : public sofa::core::topology::TopologicalMapping
{
public:
    SOFA_CLASS(AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping,sofa::core::topology::TopologicalMapping);
	typedef std::pair<size_t,SReal>  WeightedDOF;
	typedef sofa::helper::vector<std::pair<size_t,SReal> > WeightedDOFArray;
	typedef sofa::helper::vector<WeightedDOFArray> SeqWeightedDOFArray;
	typedef sofa::core::topology::Topology::Tetrahedron Tetrahedron;
	typedef sofa::defaulttype::Vec<4,SReal> Vec4;
	typedef sofa::core::topology::Topology::TetraID TetraID;
	typedef sofa::core::topology::Topology::Tetra Tetra;
	typedef sofa::core::topology::Topology::Point Point;
	typedef sofa::core::topology::Topology::Triangle Triangle;
	typedef sofa::core::topology::Topology::Edge Edge;
	typedef sofa::core::topology::Topology::Quad Quad;
	typedef sofa::core::topology::Topology::Hexahedron Hexahedron;
	typedef sofa::core::topology::BaseMeshTopology::EdgesInTriangle EdgesInTriangle;
	typedef sofa::core::topology::BaseMeshTopology::EdgesInTetrahedron EdgesInTetrahedron;
	typedef sofa::core::topology::BaseMeshTopology::EdgesInQuad EdgesInQuad;
	typedef sofa::core::topology::BaseMeshTopology::EdgesInHexahedron EdgesInHexahedron;
	typedef sofa::core::topology::BaseMeshTopology::TrianglesInTetrahedron TrianglesInTetrahedron;
	typedef sofa::core::topology::BaseMeshTopology::TetrahedraAroundTriangle TetrahedraAroundTriangle;
	template<typename D, typename E>  friend class sofa::component::mapping::AdaptiveBezierTetra2BezierTriangleMechanicalMapping;

protected:
    /** \brief Constructor.
     *
     */
    AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping();

    /** \brief Destructor.
     *
     * Does nothing.
     */
    virtual ~AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping();
public:
    /** \brief Initializes the target BaseTopology from the source BaseTopology.
     */
    virtual void init();


    /** \brief Translates the TopologyChange objects from the source to the target.
     *
     * Translates each of the TopologyChange objects waiting in the source list so that they have a meaning and
     * reflect the effects of the first topology changes on the second topology.
     *
     */
    virtual void updateTopologicalMappingTopDown();
	/// update data structure when the degree of a set of edges has changed
	void edgeChanged(const std::vector<size_t> &edges,const std::vector<size_t> &edgesInTetra);
		/// update data structure when the degree of a set of edges has changed
	void triangleChanged(const std::vector<size_t> &triangles,const std::vector<size_t> &trianglesInTetra);
    virtual unsigned int getFromIndex(unsigned int ind);

	Data<bool> flipNormals;
	Data<bool> d_useSurfaceExtrapolation;

protected:
	
	// an array providing an equivalence between each dof of the bezier triangle and the weighted dof of the bezier tetrahedron
	SeqWeightedDOFArray tetraWeightedDOFArray;
	// provide an equivalence between bezier trian vertices and bezier tetra vertices
	std::map<size_t,size_t> vertexTetra2TrianMap;
	// provide an equivalence between bezier trian vertices and bezier tetra vertices
	std::map<size_t,size_t> edgeTetra2TrianMap;
	// provide an equivalence between bezier trian vertices and bezier tetra vertices
	std::map<size_t,size_t> tetra2TrianMap;
	// for each edge indicates if the edge has been swapped from the associated edge in the tetrahedral mesh
	std::vector<bool> edgeSwapFromTetraEdge;
	// for volumetric position extrapolation store the tetrahedron associated with a given triangle
	std::map<size_t,Tetrahedron> surfaceTriangleIndexToTetrahedronMap;
	// for volumetric position extrapolation store the 2 tetrahedra associated with a given edge
	std::map<size_t, std::pair<Tetrahedron,Tetrahedron> > surfaceEdgeIndexToTetrahedronPairMap;
	// for volumetric position extrapolation store the barycentric coordinates assocated with a given edge DOF
	std::vector<std::vector<Vec4> > controlPointsBarycentricCoord;

	AdaptiveHighOrderTetrahedronSetTopologyContainer *from_btstc;
	HighOrderTriangleSetTopologyContainer *to_btstc;
	// the number of points in the output Bezier triangulation
	size_t nbDofs;

};

} // namespace topology

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_TOPOLOGY_ADAPTIVEBEZIERTETRA2BEZIERTRIANGLETOPOLOGICALMAPPING_H
