
#include "initHighOrderTopology.h"
#include "HighOrderTetra2HighOrderTriangleTopologicalMapping.h"
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>
#include "HighOrderTriangleSetTopologyContainer.h"
#include <SofaBaseTopology/TriangleSetTopologyModifier.h>
#include "HighOrderTetrahedronSetTopologyContainer.h"
#include <sofa/core/topology/TopologyChange.h>

#include <sofa/defaulttype/Vec.h>
#include <map>
#include <sofa/defaulttype/VecTypes.h>

namespace sofa
{

namespace component
{

namespace topology
{

using namespace sofa::defaulttype;

using namespace sofa::component::topology;
using namespace sofa::core::topology;

SOFA_DECL_CLASS(HighOrderTetra2HighOrderTriangleTopologicalMapping)

// Register in the Factory
int HighOrderTetra2HighOrderTriangleTopologicalMappingClass = core::RegisterObject("Special case of mapping where BezierTetrahedronSetTopology is converted to BezierTriangleSetTopology")
        .add< HighOrderTetra2HighOrderTriangleTopologicalMapping >()

        ;

// Implementation

HighOrderTetra2HighOrderTriangleTopologicalMapping::HighOrderTetra2HighOrderTriangleTopologicalMapping()
    : flipNormals(initData(&flipNormals, bool(false), "flipNormals", "Flip Normal ? (Inverse point order when creating triangle)"))
  {
}

HighOrderTetra2HighOrderTriangleTopologicalMapping::~HighOrderTetra2HighOrderTriangleTopologicalMapping()
{
}

void HighOrderTetra2HighOrderTriangleTopologicalMapping::init()
{
    //sout << "INFO_print : init HighOrderTetra2HighOrderTriangleTopologicalMapping" << sendl;

    // INITIALISATION of Bezier TRIANGULAR mesh from Bezier TETRAHEDRAL mesh :


    if (fromModel)
    {
		HighOrderTetrahedronSetTopologyContainer *from_btstc;
		fromModel->getContext()->get(from_btstc);
		if (!from_btstc) {
			serr << "Could not find an input BezierTetrahedronSetTopologyContainer"<<sendl;
		}



        if (toModel)
        {

//            sout << "INFO_print : HighOrderTetra2HighOrderTriangleTopologicalMapping - to = triangle" << sendl;

            HighOrderTriangleSetTopologyContainer *to_btstc;
            toModel->getContext()->get(to_btstc);

			if (!to_btstc) {
				serr << "Could not find an output  HighOrderTriangleSetTopologyContainer " <<sendl;
			}

            to_btstc->clear();



           TriangleSetTopologyModifier *to_tstm;
            toModel->getContext()->get(to_tstm);

            const sofa::helper::vector<Triangle> &triangleArray=fromModel->getTriangles();
            const bool flipN = flipNormals.getValue();

			// set the degree of Bezier triangle equal to that of Bezier tetra
			to_btstc->d_degree.setValue(from_btstc->getDegree());
			to_btstc->init();

			// set the number of points of Bezier triangle = number of points of Bezier tetra
            toModel->setNbPoints(from_btstc->getNbPoints());
			// initialize table of equivalence
			sofa::helper::vector <unsigned int>& Loc2GlobVec = *(Loc2GlobDataVec.beginEdit());

			Loc2GlobVec.clear();
			Glob2LocMap.clear();
			size_t rankTriangle=0;
			// set to count the number of vertices 
			std::set<size_t> vertexSet,edgeSet;
			// set the boolean indicating if the triangulation is rational
			helper::WriteOnlyAccessor<Data <HighOrderTriangleSetTopologyContainer::SeqBools> >  isRationalSpline=to_btstc->d_isRationalSpline;
            std::set<Triangle> triangleSet;
            size_t degree = from_btstc->getDegree();
            size_t j,k, globalIndex,offset;
            bool ret;
			for (unsigned int i=0; i<triangleArray.size(); ++i)
			{
				/// find triangles on the border of the tetrahedral mesh 
				TetrahedraAroundTriangle tat=fromModel->getTetrahedraAroundTriangle(i);
				if (tat.size()==1)
				{
					// add macro triangle
					Triangle t = triangleArray[i];
					if(flipN)
					{
						unsigned int tmp = t[2];
						t[2] = t[1];
						t[1] = tmp;
					}
                    triangleSet.insert(t);
					to_tstm->addTriangleProcess(t);
					// add vertices in set
					vertexSet.insert(t[0]);vertexSet.insert(t[1]);vertexSet.insert(t[2]);
                    // update local maps of control points
                    for (j = 0; j < 3; ++j) {
                        HighOrderTriangleSetTopologyContainer::ControlPointLocation cpl(t[j], std::make_pair(HighOrderTriangleSetTopologyContainer::POINT, 0));
                        to_btstc->locationToGlobalIndexMap.insert(make_pair(cpl, t[j]));
                        to_btstc->globalIndexToLocationMap.insert(make_pair(t[j], cpl));
                    }
					//  if the adjacent tetrahedron is rational then the triangle is also rational
					const bool irs=from_btstc->isRationalSpline(tat[0]);
					isRationalSpline.push_back(irs);
                    // update the 2 lookup table to get tetrahedron id from triangle id and conversely
					Loc2GlobVec.push_back(i);
					Glob2LocMap[i]=Loc2GlobVec.size()-1;
                    // update edge information
                    if (degree > 1) {
                        core::topology::BaseMeshTopology::EdgesInTriangle eit = from_btstc->getEdgesInTriangle(i);

                        for (j = 0; j < 3; ++j) {
                            if (edgeSet.count(eit[j]) == 0) {
                                Edge e = from_btstc->getEdge(eit[j]);
                                size_t v0 = std::min(e[0], e[1]);
                                size_t v1 = std::max(e[0], e[1]);
                                int realEdgeIndex = to_btstc->getEdgeIndex(v0, v1);
                                assert(realEdgeIndex >= 0);
                                Edge e2 = to_btstc->getEdge((unsigned int)realEdgeIndex);
                                for (k = 0; k < (degree - 1); ++k) {
                                   
                                    if (e[0] == e2[0]) {
                                        offset =k;
                                    }
                                    else {
                                        assert(e[0] == e2[1]);
                                        offset = degree -2 -k;
                                    }
                                    ret=from_btstc->getGlobalIndexFromLocation(HighOrderTetrahedronSetTopologyContainer::EDGE, eit[j], k, globalIndex);
                                    assert(ret);
                                    HighOrderTriangleSetTopologyContainer::ControlPointLocation cpl((size_t)realEdgeIndex, std::make_pair(HighOrderTriangleSetTopologyContainer::EDGE, offset));
                                    to_btstc->locationToGlobalIndexMap.insert(make_pair(cpl, globalIndex));
                                    to_btstc->globalIndexToLocationMap.insert(make_pair(globalIndex, cpl));
                                }
                                edgeSet.insert(eit[j]);

                            }
                        }
                    }

                    // update triangle information
                    if (degree > 2) {
                        size_t nbTriangleControlPoints = (degree - 1)*(degree - 2) / 2;
                        for (k = 0; k < nbTriangleControlPoints; ++k) {
                            ret=from_btstc->getGlobalIndexFromLocation(HighOrderTetrahedronSetTopologyContainer::TRIANGLE, i, k, globalIndex);
                            assert(ret);
                            HighOrderTriangleSetTopologyContainer::ControlPointLocation cpl(triangleSet.size()-1, std::make_pair(HighOrderTriangleSetTopologyContainer::TRIANGLE, k));
                            to_btstc->locationToGlobalIndexMap.insert(make_pair(cpl, globalIndex));
                            to_btstc->globalIndexToLocationMap.insert(make_pair(globalIndex, cpl));
                        }
                    }
                   
					rankTriangle++;
				}
			}
			// copy the weights 
			const HighOrderTetrahedronSetTopologyContainer::SeqWeights &swFrom=from_btstc->getWeightArray();

			HighOrderTriangleSetTopologyContainer::SeqWeights &wa=*(to_btstc->d_weightArray.beginEdit());
			wa.resize(swFrom.size());
			std::copy(swFrom.begin(),swFrom.end(),wa.begin());
			to_btstc->d_weightArray.endEdit();
			
			to_btstc->d_numberOfTriangularPoints.setValue(vertexSet.size());
			to_btstc->triangleDOFArray.clear(); 
			VecPointID indexArray;
			size_t i;
			for (i=0;i<(size_t)to_btstc->getNbTriangles();++i) {
				indexArray.clear();
				to_btstc->getGlobalIndexArrayOfControlPointsInTriangle(i,indexArray);
				to_btstc->triangleDOFArray.push_back(indexArray);
			}
			to_btstc->checkTopology();
			//to_tstm->propagateTopologicalChanges();
			to_tstm->notifyEndingEvent();
			//to_tstm->propagateTopologicalChanges();
			Loc2GlobDataVec.endEdit();
		}
	}
}


unsigned int HighOrderTetra2HighOrderTriangleTopologicalMapping::getFromIndex(unsigned int ind)
{

    if(fromModel->getTetrahedraAroundTriangle(ind).size()==1)
    {
        return fromModel->getTetrahedraAroundTriangle(ind)[0];
    }
    else
    {
        return 0;
    }
}

void HighOrderTetra2HighOrderTriangleTopologicalMapping::updateTopologicalMappingTopDown()
{

    // INITIALISATION of TRIANGULAR mesh from TETRAHEDRAL mesh :
//	cerr << "updateTopologicalMappingTopDown called" << endl;

    if (fromModel)
    {

        TriangleSetTopologyModifier *to_tstm;
        toModel->getContext()->get(to_tstm);

        if (toModel)
        {

            std::list<const TopologyChange *>::const_iterator itBegin=fromModel->beginChange();
            std::list<const TopologyChange *>::const_iterator itEnd=fromModel->endChange();

            //sofa::helper::vector <unsigned int>& Loc2GlobVec = *(Loc2GlobDataVec.beginEdit());

            while( itBegin != itEnd )
            {
                TopologyChangeType changeType = (*itBegin)->getChangeType();

                switch( changeType )
                {

                case core::topology::ENDING_EVENT:
                {
                    //sout << "INFO_print : Tetra2TriangleTopologicalMapping - ENDING_EVENT" << sendl;
                    to_tstm->propagateTopologicalChanges();
                    to_tstm->notifyEndingEvent();
                    to_tstm->propagateTopologicalChanges();
                    break;
                }

                case core::topology::TRIANGLESREMOVED:
                {

                    break;
                }

                case core::topology::TRIANGLESADDED:
                {
                   
                    break;
                }

                case core::topology::TETRAHEDRAADDED:
                {
                    
                    break;
                }

                case core::topology::TETRAHEDRAREMOVED:
                {
                    

                    break;

                }

                case core::topology::EDGESADDED:
                {
                    
                    break;
                }

                case core::topology::POINTSADDED:
                {

                  
                    break;
                }

                case core::topology::POINTSREMOVED:
                {
                   

                    break;
                }

                case core::topology::POINTSRENUMBERING:
                {
                    

                    break;
                }
                default:
                    // Ignore events that are not Triangle  related.
                    break;
                };

                ++itBegin;
            }
            to_tstm->propagateTopologicalChanges();
            //Loc2GlobDataVec.endEdit();
        }
    }

    return;
}


} // namespace topology

} // namespace component

} // namespace sofa
