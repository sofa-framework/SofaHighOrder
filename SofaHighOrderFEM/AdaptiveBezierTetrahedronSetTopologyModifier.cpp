
#include "AdaptiveBezierTetrahedronSetTopologyModifier.h"
#include <sofa/core/ObjectFactory.h>


namespace sofa
{

namespace component
{

namespace topology
{

using namespace std;
using namespace sofa::defaulttype;

SOFA_DECL_CLASS(AdaptiveHighOrderTetrahedronSetTopologyModifier)
int AdaptiveHighOrderTetrahedronSetTopologyModifierClass = core::RegisterObject("Adaptive High Order Tetrahedron set topology modifier")
        .add< AdaptiveHighOrderTetrahedronSetTopologyModifier >()
        ;




void AdaptiveHighOrderTetrahedronSetTopologyModifier::init()
{

	TetrahedronSetTopologyModifier::init(); // initialize the tetrahedron array
}
void AdaptiveHighOrderTetrahedronSetTopologyModifier::reinit()
{
	TetrahedronSetTopologyModifier::reinit();
}

void  AdaptiveHighOrderTetrahedronSetTopologyModifier::changeDegree(const sofa::helper::vector< unsigned int > &edges,
		const sofa::helper::vector< unsigned int > &triangles,
		const sofa::helper::vector< unsigned int > &tetrahedra)
{
	HighOrderTetrahedronDegreeChanged *e = new HighOrderTetrahedronDegreeChanged(edges,triangles,tetrahedra);
	this->addTopologyChange(e);
}



} // namespace topology

} // namespace component

} // namespace sofa
