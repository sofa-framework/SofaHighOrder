
#ifndef SOFA_HIGHORDERTOPOLOGY_HIGHORDERTOPOLOGYINDEX_H
#define SOFA_HIGHORDERTOPOLOGY_HIGHORDERTOPOLOGYINDEX_H



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


} // namespace topology

} // namespace component

} // namespace sofa

#endif
