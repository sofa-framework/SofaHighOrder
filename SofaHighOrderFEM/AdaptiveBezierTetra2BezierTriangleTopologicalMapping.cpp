//#include "initHighOrderFEM.h"
#include "AdaptiveBezierTetra2BezierTriangleTopologicalMapping.h"
#include "AdaptiveBezierTetrahedronSetTopologyModifier.h"
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>
#include "HighOrderTriangleSetTopologyContainer.h"
#include <SofaBaseTopology/TriangleSetTopologyModifier.h>
#include "AdaptiveBezierTetrahedronSetTopologyContainer.h"
#include <sofa/core/topology/TopologyChange.h>

#include <sofa/defaulttype/Vec.h>
#include <map>
#include <sofa/defaulttype/VecTypes.h>
//#include <tuple>
#include "boost/tuple/tuple.hpp"

namespace sofa
{

namespace component
{

namespace topology
{

using namespace sofa::defaulttype;

using namespace sofa::component::topology;
using namespace sofa::core::topology;

SOFA_DECL_CLASS(AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping)

// Register in the Factory
int AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMappingClass = core::RegisterObject("Special case of mapping where AdaptiveBezierTetrahedronSetTopology is converted to BezierTriangleSetTopology")
        .add< AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping >()

        ;

// Implementation
typedef AdaptiveHighOrderTetrahedronSetTopologyContainer::ControlPointLocation ControlPointLocation;
typedef std::map<AdaptiveHighOrderTetrahedronSetTopologyContainer::ControlPointLocation,size_t>::iterator Loc2GlobMapIterator;

AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping::AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping()
    : flipNormals(initData(&flipNormals, bool(false), "flipNormals", "Flip Normal ? (Inverse point order when creating triangle)"))
	, d_useSurfaceExtrapolation(initData(&d_useSurfaceExtrapolation, true,"useSurfaceExtrapolation", "indicates if surface or volumetric extrapolation is use to update the surface triangles"))
  {
}

AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping::~AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping()
{
}
void AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping::triangleChanged(const std::vector<size_t> &triangles,
																		   const std::vector<size_t> &trianglesInTetra)

{
	helper::ReadAccessor<Data <AdaptiveHighOrderTetrahedronSetTopologyContainer::SeqBezierDegree> >  triangleDegree=from_btstc->d_triangleDegree;

	size_t i,j,k,l;
	unsigned int degree=from_btstc->getDegree();
	WeightedDOF wdf;
	WeightedDOFArray wdfa(1);
	Loc2GlobMapIterator itl2g;
	ControlPointLocation cpl;
	bool useSurfaceExtrapolation=d_useSurfaceExtrapolation.getValue();


	for (i=0;i<triangles.size();++i) {
		if (triangleDegree[trianglesInTetra[i]]==degree) 
			wdfa.resize(1);
		else 
			wdfa.resize(3);
		size_t triangleRank=to_btstc->getNumberOfTriangularPoints()+edgeTetra2TrianMap.size()*(degree-1)+
			triangles[i]*(degree-1)*(degree-2)/2;
		Triangle torg =from_btstc->getTriangle(trianglesInTetra[i]);


		for (l=0,j=1;j<(size_t)(degree-1);++j) {
			for (k=1;k<(degree-j);++k,++l,++triangleRank) {


				if (triangleDegree[trianglesInTetra[i]]==degree)  {
					cpl=ControlPointLocation(trianglesInTetra[i],std::make_pair(HighOrderTetrahedronSetTopologyContainer::TRIANGLE,l));
					// get the index of the associated DOF from the map.
					itl2g=from_btstc->locationToGlobalIndexMap.find(cpl);
					// check that the DOF exists.
					assert(itl2g!=from_btstc->locationToGlobalIndexMap.end());
					// set the mapped dof
					wdf=WeightedDOF((*itl2g).second,(SReal)1.0f);
					wdfa[0]=wdf;
					vertexTetra2TrianMap.insert(std::make_pair((*itl2g).second,triangleRank));
					tetraWeightedDOFArray[triangleRank]=wdfa;
				} else {
					// update map
					WeightedDOFArray tmp=tetraWeightedDOFArray[triangleRank];
					assert(tmp.size()==1);
					std::map<size_t,size_t>::iterator itv=vertexTetra2TrianMap.find(tmp[0].first);
//					if(itv!=vertexTetra2TrianMap.end())
//						vertexTetra2TrianMap.erase(itv);
					// the DOF are linear mapped from the 3 extremities of the triangle
					wdf=WeightedDOF(torg[0],(SReal)(j)/degree);
					wdfa[0]=wdf;
					wdf=WeightedDOF(torg[1],(SReal)(k)/degree);
					wdfa[1]=wdf;
					wdf=WeightedDOF(torg[2],(SReal)(degree-j-k)/degree);
					wdfa[2]=wdf;
					tetraWeightedDOFArray[triangleRank]=wdfa;
					if (useSurfaceExtrapolation==false) {
						// update the extrapolation weights
						// get the barycentric coordinates from the tetrahedron and update edgeControlPointsBarycentricCoord
						std::map<std::pair<size_t,size_t>,Vec4 > ::iterator ittr=
							from_btstc->triangleControlPointsBarycentricCoord.find(std::make_pair(trianglesInTetra[i],l));
						std::pair<size_t,size_t> ind2(triangles[i],l);
						std::vector<Vec4> tbcArray(1);
						tbcArray[0]=(*ittr).second;
						controlPointsBarycentricCoord[triangleRank]=tbcArray;
					}
				}
			}
		}
	}
}
void AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping::edgeChanged(const std::vector<size_t> &edges,
																	   const std::vector<size_t> &edgesInTetra
																	   )
{
	helper::ReadAccessor<Data <AdaptiveHighOrderTetrahedronSetTopologyContainer::SeqBezierDegree> >  edgeDegree=from_btstc->d_edgeDegree;

	size_t i,k;
	unsigned int degree=from_btstc->getDegree();
	WeightedDOF wdf;
	WeightedDOFArray wdfa(1);
	Loc2GlobMapIterator itl2g;
	bool useSurfaceExtrapolation=d_useSurfaceExtrapolation.getValue();
	
	for (i=0;i<edges.size();++i) {
		if (edgeDegree[edgesInTetra[i]]==degree) 
			wdfa.resize(1);
		else 
			wdfa.resize(2);
		size_t edgeRank=to_btstc->getNumberOfTriangularPoints()+edges[i]*(degree-1);
		if (edgeDegree[edgesInTetra[i]]==degree)  {
			// edge is full degree therefore there is 1 tetra DOF for each triangle DOF
			for (k=1;k<(size_t)(degree);++k,++edgeRank) {
				ControlPointLocation cpl;
				// test orientation of edge in the tetrahedron
				if (edgeSwapFromTetraEdge[edges[i]]==false)
					cpl=ControlPointLocation(edgesInTetra[i],std::make_pair(HighOrderTetrahedronSetTopologyContainer::EDGE,k-1));
				else {
					cpl=ControlPointLocation(edgesInTetra[i],std::make_pair(HighOrderTetrahedronSetTopologyContainer::EDGE,degree-k-1));
				}
				// get the index of the associated DOF from the map.
				itl2g=from_btstc->locationToGlobalIndexMap.find(cpl);
				// check that the DOF exists.
				assert(itl2g!=from_btstc->locationToGlobalIndexMap.end());

				
				// set the mapped dof
				wdf=WeightedDOF((*itl2g).second,(SReal)1.0f);
				wdfa[0]=wdf;
				vertexTetra2TrianMap.insert(std::make_pair((*itl2g).second,edgeRank));
				tetraWeightedDOFArray[edgeRank]=wdfa;
			} 
		} else {
			// the DOF are linear mapped from the 2 extremities of the edge
			Edge e=from_btstc->getEdge(edgesInTetra[i]);
			if (edgeSwapFromTetraEdge[edges[i]]) {
				unsigned int val=e[0];
				e[0]=e[1];
				e[1]=val;
			}
			for (k=1;k<(size_t)(degree);++k,++edgeRank) {
				// update map
				WeightedDOFArray tmp=tetraWeightedDOFArray[edgeRank];
				assert(tmp.size()==1);
				std::map<size_t,size_t>::iterator itv=vertexTetra2TrianMap.find(tmp[0].first);
			/*	if (itv!=vertexTetra2TrianMap.end()) {
					std::cerr<<"erasing"<<std::endl;
					vertexTetra2TrianMap.erase(itv);
				} */
				
				
				wdf=WeightedDOF(e[0],(SReal)(degree-k)/degree);
				wdfa[0]=wdf;
				wdf=WeightedDOF(e[1],(SReal)(k)/degree);
				wdfa[1]=wdf;
				tetraWeightedDOFArray[edgeRank]=wdfa;

				// update the extrapolation weights
				size_t ind=k-1;
				if (edgeSwapFromTetraEdge[edges[i]]) {
					ind=degree-k-1;
				}
				if (useSurfaceExtrapolation==false) {
					// get the barycentric coordinates from the tetrahedron and update edgeControlPointsBarycentricCoord
					std::map<std::pair<size_t,size_t>,std::pair<Vec4,Vec4> > ::iterator ite=
						from_btstc->edgeControlPointsBarycentricCoord.find(std::make_pair(edgesInTetra[i],ind));
					std::pair<size_t,size_t> ind2(edges[i],k-1);
					std::vector<Vec4> tbcArray(2);
					tbcArray[0]=(*ite).second.first;
					tbcArray[1]=(*ite).second.second;
					controlPointsBarycentricCoord[edgeRank]=tbcArray;
				}
			}
		}

	}

}
void AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping::init()
{
    //sout << "INFO_print : init BezierTetra2BezierTriangleTopologicalMapping" << sendl;

    // INITIALISATION of Bezier TRIANGULAR mesh from Bezier TETRAHEDRAL mesh :


    if (fromModel)
    {
		
		fromModel->getContext()->get(from_btstc);
		if (!from_btstc) {
			serr << "Could not find an input HighOrderTetrahedronSetTopologyContainer"<<sendl;
		}



        if (toModel)
        {

//            sout << "INFO_print : BezierTetra2BezierTriangleTopologicalMapping - to = triangle" << sendl;

            
            toModel->getContext()->get(to_btstc);

			if (!to_btstc) {
				serr << "Could not find an output  HighOrderTriangleSetTopologyContainer " <<sendl;
			}

			to_btstc->clear();
			// set the degree of Bezier triangle equal to that of Bezier tetra
			to_btstc->d_degree.setValue(from_btstc->getDegree());


           TriangleSetTopologyModifier *to_tstm;
            toModel->getContext()->get(to_tstm);

            const sofa::helper::vector<Triangle> &triangleArray=fromModel->getTriangles();
            const bool flipN = flipNormals.getValue();

			// initialize table of equivalence
			helper::WriteOnlyAccessor<Data <sofa::helper::vector <unsigned int> > > Loc2GlobVec = Loc2GlobDataVec; 

			Loc2GlobVec.clear();
			Glob2LocMap.clear();
			size_t rankTriangle=0;
			// set to count the number of vertices 
			std::set<size_t> vertexSet;
			// count the number of vertex and store equivalence map.
			std::set<size_t>::iterator ite;
			// parse the edges on surface to get the vertices on surface
			for (ite=from_btstc->edgeOnSurfaceSet.begin();ite!=from_btstc->edgeOnSurfaceSet.end();++ite) {
				Edge e=from_btstc->getEdge(*ite);
				vertexSet.insert(e[0]);
				vertexSet.insert(e[1]);
			}
			//  weights from input tetrahedron and output bezier triangulation
			const HighOrderTetrahedronSetTopologyContainer::SeqWeights &swFrom=from_btstc->getWeightArray();
			helper::WriteOnlyAccessor<Data <HighOrderTriangleSetTopologyContainer::SeqWeights> >  swTo=to_btstc->d_weightArray;

			// now parse vertex on surface
			std::set<size_t>::iterator itv;
			WeightedDOF wdf;
			WeightedDOFArray wdfa(1);
			// fill vertexTrian2TetraMap and tetraWeightedDOFArray for vertices
			tetraWeightedDOFArray.clear();
			vertexTetra2TrianMap.clear();
			for (itv=vertexSet.begin();itv!=vertexSet.end();++itv) {
				wdf=WeightedDOF(*itv,(SReal)1.0f);
				wdfa[0]=wdf;
				// update map to get the index of triangle vertex from tetrahedrization vertex
				vertexTetra2TrianMap[*itv]=tetraWeightedDOFArray.size();
				// add DOF description
				tetraWeightedDOFArray.push_back(wdfa);
				// add weight
				swTo.push_back(swFrom[*itv]);
			}
			// set the number of triangular points = vertices of the triangulation
			to_btstc->d_numberOfTriangularPoints.setValue(vertexSet.size());
		
			// set the boolean indicating if the triangulation is rational
			helper::WriteOnlyAccessor<Data <HighOrderTriangleSetTopologyContainer::SeqBools> >  isRationalSpline=to_btstc->d_isRationalSpline;
			std::map<size_t,size_t>::iterator itvtet; 

			std::vector<Edge> edgeArray;
			std::map<Edge,size_t> edgeMap;
			std::map<Edge,size_t>::iterator item;
			std::vector<Triangle> newTriangleArray;
			unsigned int i,j,k;
			unsigned int degree=from_btstc->getDegree();

			nbDofs=vertexSet.size()+(degree-1)*from_btstc->edgeOnSurfaceSet.size()+
				(degree-1)*(degree-2)*from_btstc->triangleOnSurfaceSet.size()/2;
			to_btstc->setNbPoints(nbDofs);
			controlPointsBarycentricCoord.resize(nbDofs);

			helper::ReadAccessor<Data <AdaptiveHighOrderTetrahedronSetTopologyContainer::SeqBezierDegree> >  edgeDegree=from_btstc->d_edgeDegree;
			helper::ReadAccessor<Data <AdaptiveHighOrderTetrahedronSetTopologyContainer::SeqBezierDegree> >  triangleDegree=from_btstc->d_triangleDegree;
			Loc2GlobMapIterator itl2g;

			for (i=0; i<triangleArray.size(); ++i)
			{
				/// find triangles on the border of the tetrahedral mesh 
				TetrahedraAroundTriangle tat=fromModel->getTetrahedraAroundTriangle(i);
				if (tat.size()==1)
				{
					// store the tetrahedron index associated with each triangle in a  map
					tetra2TrianMap[tat[0]]=Loc2GlobVec.size();
					// add macro triangle torg = triangle in bezier tetrahedra t = triangle in bezier triangle
					Triangle torg = triangleArray[i];
					Triangle t=torg;

					if(flipN)
					{
						unsigned int tmp = t[2];
						t[2] = t[1];
						t[1] = tmp;
					}
					// renumber triangle with new triangle vertices
					for (j=0;j<3;++j) {
						itvtet=vertexTetra2TrianMap.find(t[j]);
						assert (itvtet!=vertexTetra2TrianMap.end());
						t[j]=(*itvtet).second;
					}
					// add triangle
					to_tstm->addTriangleProcess(t);
					newTriangleArray.push_back(t);
					
					// equivalence between triangles indices between Bezier tetrahedra and Bezier triangles
					Loc2GlobVec.push_back(i);
					Glob2LocMap[i]=Loc2GlobVec.size()-1;
					//  if the adjacent tetrahedron is rational then the triangle is also rational
					const bool irs=from_btstc->isRationalSpline(tat[0]);
					isRationalSpline.push_back(irs);
					// now add description of edge control point in Bezier triangulation
					if (degree>1) {

						EdgesInTriangle eit=from_btstc->getEdgesInTriangle(i);
						// now handles edge if necessary
						for (j=0;j<3;++j) {
							Edge eOrg=from_btstc->getEdge(eit[j]);
							Edge e(t[(j+1)%3],t[(j+2)%3]);
							Edge se=Edge(std::min(e[0],e[1]),std::max(e[0],e[1]));
							if ((item=edgeMap.find(se))==edgeMap.end()) {
								// new edge : update map and array
								edgeMap.insert(std::make_pair(se,edgeArray.size()));
								edgeArray.push_back(e);
								// update edge map with tetra edges
								edgeTetra2TrianMap[eit[j]]=edgeArray.size()-1;
								// set the weight of the DOF as the one from the 
								// test degree of edge
								if (edgeDegree[eit[j]]==degree) 
									wdfa.resize(1);
								else 
									wdfa.resize(2);
								// store a boolean indicating if the Bezier triangulation edge
								// has the same orientation  or not as the Bezier tetrahedral edge
								if (eOrg[0]==torg[(j+1)%3]) {
									edgeSwapFromTetraEdge.push_back(false);
								} else
									edgeSwapFromTetraEdge.push_back(true);
								// update surfaceEdgeIndexToTetrahedronPairMap
								std::map<size_t, std::pair<size_t,size_t> >::iterator ite=
									from_btstc->edgeOnSurfaceMap.find(eit[j]);
								assert(ite!=from_btstc->edgeOnSurfaceMap.end());
								Tetrahedron tet1=from_btstc->getTetrahedron((*ite).second.first);
								Tetrahedron tet2=from_btstc->getTetrahedron((*ite).second.second);
								size_t ind=edgeArray.size()-1;
								surfaceEdgeIndexToTetrahedronPairMap[ind]=std::make_pair(tet1,tet2);
								for (k=1;k<(size_t)(degree);++k) {
									ControlPointLocation cpl;
									// test orientation of edge in the tetrahedron
									if (eOrg[0]==torg[(j+1)%3]) {
										cpl=ControlPointLocation(eit[j],std::make_pair(HighOrderTetrahedronSetTopologyContainer::EDGE,k-1));
									} else {
										cpl=ControlPointLocation(eit[j],std::make_pair(HighOrderTetrahedronSetTopologyContainer::EDGE,degree-k-1));
									}
									// get the index of the associated DOF from the map.
									itl2g=from_btstc->locationToGlobalIndexMap.find(cpl);
									// check that the DOF exists.
									assert(itl2g!=from_btstc->locationToGlobalIndexMap.end());
									// fill edgeControlPointsBarycentricCoord with dummy value
									std::pair<size_t,size_t> ind2(ind,k-1);
	
									// specify weight
									swTo.push_back(swFrom[(*itl2g).second]);
									if (edgeDegree[eit[j]]==degree)  {
										// set the mapped dof
										wdf=WeightedDOF((*itl2g).second,(SReal)1.0f);
										wdfa[0]=wdf;
										vertexTetra2TrianMap.insert(std::make_pair((*itl2g).second,tetraWeightedDOFArray.size()));
										tetraWeightedDOFArray.push_back(wdfa);

									} else {
										// the DOF are linear mapped from the 2 extremities of the edge
										wdf=WeightedDOF(torg[(j+1)%3],(SReal)(degree-k)/degree);
										wdfa[0]=wdf;
										wdf=WeightedDOF(torg[(j+2)%3],(SReal)(k)/degree);
										wdfa[1]=wdf;
										tetraWeightedDOFArray.push_back(wdfa);
									}

									
								}
							}
						}
					}
				}
			}
			if (degree >2) {
				size_t l;
				for (i=0; i<triangleArray.size(); ++i)
				{
					/// find triangles on the border of the tetrahedral mesh 
					TetrahedraAroundTriangle tat=fromModel->getTetrahedraAroundTriangle(i);
					if (tat.size()==1)
					{
						Triangle torg = triangleArray[i];
						ControlPointLocation cpl;
						if (triangleDegree[i]==degree)  
							wdfa.resize(1);
						else
							wdfa.resize(3);
						// update edgeOnSurfaceMap
						std::map<size_t, size_t>::iterator ittr=
							from_btstc->trianglesOnSurfaceMap.find(i);
						assert(ittr!=from_btstc->trianglesOnSurfaceMap.end());
						Tetrahedron tet=from_btstc->getTetrahedron((*ittr).second);
						size_t ind=Glob2LocMap[i];
						surfaceTriangleIndexToTetrahedronMap[ind]=tet;
						for (l=0,j=1;j<(size_t)(degree-1);++j) {
							for (k=1;k<(degree-j);++k,++l) {
								cpl=ControlPointLocation(i,std::make_pair(HighOrderTetrahedronSetTopologyContainer::TRIANGLE,l));
								// get the index of the associated DOF from the map.
								itl2g=from_btstc->locationToGlobalIndexMap.find(cpl);
								// check that the DOF exists.
								assert(itl2g!=from_btstc->locationToGlobalIndexMap.end());
								// specify weight
								swTo.push_back(swFrom[(*itl2g).second]);
								// fill edgeControlPointsBarycentricCoord with dummy value
								std::pair<size_t,size_t> ind2(ind,l);
								if (triangleDegree[i]==degree)  {
									// set the mapped dof
									wdf=WeightedDOF((*itl2g).second,(SReal)1.0f);
									wdfa[0]=wdf;
									vertexTetra2TrianMap.insert(std::make_pair((*itl2g).second,tetraWeightedDOFArray.size()));
									tetraWeightedDOFArray.push_back(wdfa);
								} else {
									// the DOF are linear mapped from the 3 extremities of the triangle
									wdf=WeightedDOF(torg[0],(SReal)(j)/degree);
									wdfa[0]=wdf;
									wdf=WeightedDOF(torg[1],(SReal)(k)/degree);
									wdfa[1]=wdf;
									wdf=WeightedDOF(torg[2],(SReal)(degree-j-k)/degree);
									wdfa[2]=wdf;
									tetraWeightedDOFArray.push_back(wdfa);
								}

							}
						}
					}


				}
			}

			to_btstc->checkTopology();
			//to_tstm->propagateTopologicalChanges();
			to_tstm->notifyEndingEvent();
			//to_tstm->propagateTopologicalChanges();

			to_btstc->init();
			
		}
	}


}


unsigned int AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping::getFromIndex(unsigned int ind)
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

void AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping::updateTopologicalMappingTopDown()
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
				case core::topology::TOPOLOGYCHANGE_LASTID:
                {
					const HighOrderTetrahedronDegreeChanged *btdc=dynamic_cast<const HighOrderTetrahedronDegreeChanged *>(*itBegin);
					assert(btdc!=NULL);
					if (btdc) {
						std::map<size_t,size_t>::iterator itchanged;
						const sofa::helper::vector<unsigned int> &edgeModified=btdc->getEdgeArray();
						std::vector<size_t> edges,edgesInTetra;
						size_t i;
						for (i=0;i<edgeModified.size();++i) {
							itchanged=edgeTetra2TrianMap.find(edgeModified[i]);
							
							if (itchanged!=edgeTetra2TrianMap.end()) {
								// edge changed on Bezier triangle
								edges.push_back((*itchanged).second);
								edgesInTetra.push_back((*itchanged).first);
							}
						}
						edgeChanged(edges,edgesInTetra);
						const sofa::helper::vector<unsigned int> &triangleModified=btdc->getTriangleArray();
						std::vector<size_t> triangles,trianglesInTetra;
						std::map<unsigned int,unsigned int>::iterator itTrianChanged;
						for (i=0;i<triangleModified.size();++i) {
							itTrianChanged=Glob2LocMap.find(triangleModified[i]);
							
							if (itTrianChanged!=Glob2LocMap.end()) {
								// edge changed on Bezier triangle
								triangles.push_back((*itTrianChanged).second);
								trianglesInTetra.push_back((*itTrianChanged).first);
							}
						}
						triangleChanged(triangles,trianglesInTetra);
						// check consistency
#ifndef NDEBUG
						std::map<size_t,size_t>::iterator itv=vertexTetra2TrianMap.begin();
						for (;itv!=vertexTetra2TrianMap.end();++itv) {
							assert( tetraWeightedDOFArray[(*itv).second].size()==1);
							assert( tetraWeightedDOFArray[(*itv).second][0].first==(*itv).first);
						}
						for (i=0;i<tetraWeightedDOFArray.size();++i) {
							if (tetraWeightedDOFArray[i].size()==1) {
								assert(vertexTetra2TrianMap.find(tetraWeightedDOFArray[i][0].first)!=vertexTetra2TrianMap.end());
								itv=vertexTetra2TrianMap.find(tetraWeightedDOFArray[i][0].first);
								assert((*itv).second==i);
							}
						}
#endif
						to_tstm->propagateTopologicalChanges();
					}
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
						const PointsRemoved *pr=dynamic_cast<const PointsRemoved *>(*itBegin);
						if (pr) {
							 const sofa::helper::vector<unsigned int> removedPointsInTetra=pr->getArray();
							 sofa::helper::vector<unsigned int> updatedPointsInTrian;
							 sofa::helper::vector< sofa::helper::vector< unsigned int > > ancestorsList;
							 sofa::helper::vector< sofa::helper::vector< double > > baryCoefsList;
							 size_t i;
							 // parse points to only keep those on the surface
							 std::map<size_t,size_t>::iterator itv=vertexTetra2TrianMap.begin();
							 std::vector<boost::tuple<size_t,size_t,size_t> > swapArray;
							 size_t trianIndex;
							
							 sofa::helper::vector< unsigned int > al(1);
							 sofa::helper::vector< double > bcl(1);
							 bcl[0]=1.0;
							 for (i=0;i<removedPointsInTetra.size();++i) {
								 itv=vertexTetra2TrianMap.find(removedPointsInTetra[i]);
								 if (itv!=vertexTetra2TrianMap.end()) {
									 // this point is on the surface
									 updatedPointsInTrian.push_back((*itv).second);

									 al[0]=(*itv).first;
									 ancestorsList.push_back(al);
									 baryCoefsList.push_back(bcl);
								 }
							 }
							 // now send a topological event to component indicating that points have  moved
							 to_tstm->PointSetTopologyModifier::movePointsProcess(updatedPointsInTrian,ancestorsList,baryCoefsList,false);
							 to_tstm->propagateTopologicalChanges();
							 // now take into account the renumbering of tetra vertices 
							 size_t lastPoint=from_btstc->getNbPoints()-1;
							 for (i=0;i<removedPointsInTetra.size();++i,--lastPoint) {
								 itv=vertexTetra2TrianMap.find(removedPointsInTetra[i]);
								 if (itv!=vertexTetra2TrianMap.end()) {
									 trianIndex=(*itv).second;
									 vertexTetra2TrianMap.erase(itv);
						//			 vertexTetra2TrianMap[lastPoint]=trianIndex;
									 // update tetraWeightedDOFArray
									 assert(tetraWeightedDOFArray[trianIndex].size()==1);
									 WeightedDOF wdf;
									 WeightedDOFArray wdfa;
									 wdfa=tetraWeightedDOFArray[trianIndex];
									 assert(wdfa[0].first==removedPointsInTetra[i]);
									 wdfa[0].first=lastPoint;
									 tetraWeightedDOFArray[trianIndex]=wdfa;
								 }

								 // now process the last point
								 if (removedPointsInTetra[i]!=lastPoint) {
									 itv=vertexTetra2TrianMap.find(lastPoint);
									 if (itv!=vertexTetra2TrianMap.end()) {
										 // this point has changed index
										 trianIndex=(*itv).second;
										 vertexTetra2TrianMap.erase(itv);
										 vertexTetra2TrianMap[removedPointsInTetra[i]]=trianIndex;
										 WeightedDOF wdf;
										 WeightedDOFArray wdfa;
										 wdfa=tetraWeightedDOFArray[trianIndex];
										 assert(wdfa[0].first==lastPoint);
										 wdfa[0].first=removedPointsInTetra[i];
										 tetraWeightedDOFArray[trianIndex]=wdfa;
									 }
								 }
							 }

						}
		//				std::cerr<<" got points removed2"<<std::endl;
						break;
                }

                case core::topology::POINTSRENUMBERING:
                {
                    

                    break;
                }
                default:
                    // Ignore events that are not Triangle  related.
					 to_tstm->propagateTopologicalChanges();
                    break;
                };

                ++itBegin;
            }
        //    to_tstm->propagateTopologicalChanges();
            //Loc2GlobDataVec.endEdit();
        }
    }

    return;
}


} // namespace topology

} // namespace component

} // namespace sofa
