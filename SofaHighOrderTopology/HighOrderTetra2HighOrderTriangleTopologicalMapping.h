#ifndef SOFA_HIGHORDERTOPOLOGY_BEZIERTETRA2BEZIERTRIANGLETOPOLOGICALMAPPING_H
#define SOFA_HIGHORDERTOPOLOGY_BEZIERTETRA2BEZIERTRIANGLETOPOLOGICALMAPPING_H

#include <sofa/core/topology/TopologicalMapping.h>

#include <sofa/defaulttype/Vec.h>
#include <map>

#include <sofa/core/BaseMapping.h>

namespace sofa
{

namespace component
{

namespace topology
{


/**
 * This class, called HighOrderTetra2HighOrderTriangleTopologicalMapping, is a specific implementation of the interface TopologicalMapping where :
 *
 * INPUT TOPOLOGY = HighOrderTetrahedronSetTopology
 * OUTPUT TOPOLOGY = HighOrderTriangleSetTopology, as the boundary of the INPUT TOPOLOGY
 *
 * HighOrderTetra2HighOrderTriangleTopologicalMapping class is templated by the pair (INPUT TOPOLOGY, OUTPUT TOPOLOGY)
 *
*/

class SOFA_HIGHORDER_TOPOLOGY_API HighOrderTetra2HighOrderTriangleTopologicalMapping : public sofa::core::topology::TopologicalMapping
{
public:
    SOFA_CLASS(HighOrderTetra2HighOrderTriangleTopologicalMapping,sofa::core::topology::TopologicalMapping);
	typedef sofa::core::topology::Topology::Tetrahedron Tetrahedron;
	typedef sofa::core::topology::Topology::PointID			         PointID;
	typedef sofa::helper::vector<PointID>         VecPointID;
	typedef sofa::core::topology::Topology::TetraID TetraID;
	typedef sofa::core::topology::Topology::Tetra Tetra;
	typedef sofa::core::topology::Topology::Edge Edge;
	typedef sofa::core::topology::Topology::Triangle Triangle;
	typedef sofa::core::topology::BaseMeshTopology::SeqTriangles SeqTriangles;
	typedef sofa::core::topology::BaseMeshTopology::TrianglesInTetrahedron TrianglesInTetrahedron;
	typedef sofa::core::topology::BaseMeshTopology::TetrahedraAroundTriangle TetrahedraAroundTriangle;
protected:
    /** \brief Constructor.
     *
     */
    HighOrderTetra2HighOrderTriangleTopologicalMapping();

    /** \brief Destructor.
     *
     * Does nothing.
     */
    virtual ~HighOrderTetra2HighOrderTriangleTopologicalMapping();
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

    virtual unsigned int getFromIndex(unsigned int ind);
protected:
    Data<bool> flipNormals;


    std::vector<unsigned int> addedTriangleIndex;
};

} // namespace topology

} // namespace component

} // namespace sofa

#endif // SOFA_HIGHORDERTOPOLOGY_BEZIERTETRA2BEZIERTRIANGLETOPOLOGICALMAPPING_H
