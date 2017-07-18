
#ifndef SOFA_HIGHORDERTOPOLOGY_HIGHORDERTOPOLOGYINDEX_H
#define SOFA_HIGHORDERTOPOLOGY_HIGHORDERTOPOLOGYINDEX_H
#include <SofaBaseTopology/CommonAlgorithms.h>
#include "initHighOrderTopology.h"

namespace sofa
{

namespace component
{

namespace topology
{


typedef unsigned char HighOrderDegreeType;
typedef sofa::defaulttype::Vec<2,HighOrderDegreeType>		EdgeIndexVector;
typedef sofa::defaulttype::Vec<2,HighOrderDegreeType>		QuadrilateralIndexVector;
typedef sofa::defaulttype::Vec<3,HighOrderDegreeType>		HexahedronIndexVector;
typedef sofa::defaulttype::Vec<3,HighOrderDegreeType>		TriangleIndexVector;
typedef sofa::defaulttype::Vec<4,HighOrderDegreeType>		TetrahedronIndexVector;

 size_t  SOFA_HIGHORDER_TOPOLOGY_API factorialTVI(TriangleIndexVector tvi);

 size_t  SOFA_HIGHORDER_TOPOLOGY_API factorialTVI(TetrahedronIndexVector tvi);

} // namespace topology

} // namespace component

} // namespace sofa

#endif
