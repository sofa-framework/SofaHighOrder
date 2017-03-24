
#ifndef SOFA_HIGHORDERTOPOLOGY_MESH2BEZIERTOPOLOGICALMAPPING_H
#define SOFA_HIGHORDERTOPOLOGY_MESH2BEZIERTOPOLOGICALMAPPING_H
#include "initHighOrderTopology.h"

#include <SofaTopologyMapping/Mesh2PointTopologicalMapping.h>



namespace sofa
{
namespace component
{
namespace topology
{

/**
 * This class, called Mesh2HighOrderTopologicalMapping, is a specific implementation of the interface TopologicalMapping where :
 *
 * INPUT TOPOLOGY = any Tetrahedral or triangular MeshTopology
 * OUTPUT TOPOLOGY = A HighOrderTetrahedronSetTopology or BezierTriangleSetTopology as a tesselated version of the input mesh 
 *
 * This Topological mapping is a specific implementation of the Mesh2PointTopologicalMapping with a small overhead
 *
 * Mesh2HighOrderTopologicalMapping class is templated by the pair (INPUT TOPOLOGY, OUTPUT TOPOLOGY)
 *
*/

    class SOFA_HIGHORDER_TOPOLOGY_API Mesh2HighOrderTopologicalMapping : public sofa::component::topology::Mesh2PointTopologicalMapping
{
public:
    SOFA_CLASS(Mesh2HighOrderTopologicalMapping,sofa::component::topology::Mesh2PointTopologicalMapping);
protected:
    /** \brief Constructor.
     *
     */
    Mesh2HighOrderTopologicalMapping ();

    /** \brief Destructor.
     *
         * Does nothing.
         */
    virtual ~Mesh2HighOrderTopologicalMapping() {};
public:
    /** \brief Initializes the target BaseTopology from the source BaseTopology.
     */
    virtual void init();

    /// Fills pointBaryCoords, edgeBaryCoords, triangleBaryCoords and tetraBaryCoords so as to create a Bezier Tetrahedron mesh of a given order
	Data < unsigned int > bezierTetrahedronDegree;
	/// Fills pointBaryCoords, edgeBaryCoords, triangleBaryCoords so as to create a Bezier Triangle mesh of a given order
	Data < unsigned int > bezierTriangleDegree;
};

} // namespace topology
} // namespace component
} // namespace sofa

#endif // SOFA_HIGHORDERTOPOLOGY_MESH2BEZIERTOPOLOGICALMAPPING_H
