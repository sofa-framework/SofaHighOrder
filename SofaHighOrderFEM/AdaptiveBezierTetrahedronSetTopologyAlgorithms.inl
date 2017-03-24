
#ifndef SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERTETRAHEDRONSETTOPOLOGYALGORITHMS_INL
#define SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERTETRAHEDRONSETTOPOLOGYALGORITHMS_INL

#include "AdaptiveBezierTetrahedronSetTopologyAlgorithms.h"
#include "AdaptiveBezierTetrahedronSetTopologyContainer.h"
#include "AdaptiveBezierTetrahedronSetTopologyModifier.h"
#include <BezierTetrahedronSetGeometryAlgorithms.h>
#include <SofaBaseTopology/TetrahedronSetTopologyModifier.h>
#include <sofa/core/visual/VisualParams.h>


namespace sofa
{

namespace component
{

namespace topology
{
	typedef std::map<AdaptiveHighOrderTetrahedronSetTopologyContainer::ControlPointLocation,size_t>::iterator Loc2GlobMapIterator;
	typedef std::map<size_t,AdaptiveHighOrderTetrahedronSetTopologyContainer::ControlPointLocation>::iterator Glob2LocMapIterator;
	typedef AdaptiveHighOrderTetrahedronSetTopologyContainer::ControlPointLocation ControlPointLocation;
	typedef  AdaptiveHighOrderTetrahedronSetTopologyContainer::Vec4 Vec4;


	template<class DataTypes>
 AdaptiveBezierTetrahedronSetTopologyAlgorithms< DataTypes >::AdaptiveBezierTetrahedronSetTopologyAlgorithms()
        : TetrahedronSetTopologyAlgorithms<DataTypes>()
	{


}

template<class DataTypes>
void AdaptiveBezierTetrahedronSetTopologyAlgorithms< DataTypes >::init()
{
    TetrahedronSetTopologyAlgorithms< DataTypes >::init();
    this->getContext()->get(m_adaptiveContainer);
	this->getContext()->get(m_adaptiveModifier);
    this->getContext()->get(m_bezierGeometryAlgorithms);


	bool useSurfaceExtrapolation=m_adaptiveContainer->d_useSurfaceExtrapolation.getValue();

	object = this->getContext()->core::objectmodel::BaseContext::template get< core::behavior::MechanicalState< DataTypes > >();

	size_t i;


	// store the rest position from the tetrahedron mechanical state.
	size_t globalIndex,elementOffset,elementIndex;
	HighOrderTetrahedronPointLocation location;
    const VecCoord &restPosition = object->read(core::ConstVecCoordId::restPosition())->getValue();
	for (i=0;i<m_adaptiveContainer->getNbPoints();++i) {
		m_adaptiveContainer->getLocationFromGlobalIndex(i,location,elementIndex,elementOffset);
		ControlPointLocation cpl(elementIndex,std::make_pair(location,elementOffset));
		restPositionMap[cpl]=restPosition[i];
		offsetPositionMap[cpl]=Coord();
	}
	// fill restPositionArray
	restPositionArray.resize(m_adaptiveContainer->getNumberOfTetrahedralPoints());
	for (i=0;i<restPositionArray.size();++i) {
		restPositionArray[i]=restPosition[i];
	}

//	for (i=0;i<btstc->getNumberOfTriangularPoints();++i) {
//		restPositionArray.push_back(restPosition[topoMap->tetraWeightedDOFArray[i][0].first]);
//	}
}

/// computes the barycentric coordinates of a point with respect to a tetrahedron
template<class DataTypes>
typename AdaptiveHighOrderTetrahedronSetTopologyContainer::Vec4 getBarycentricCoordinates(
typename	AdaptiveBezierTetrahedronSetTopologyAlgorithms< DataTypes >::Coord &pos,
typename	AdaptiveBezierTetrahedronSetTopologyAlgorithms< DataTypes >::Coord tetraPos[4]) {
		AdaptiveHighOrderTetrahedronSetTopologyContainer::Vec4 v;

        typename AdaptiveBezierTetrahedronSetTopologyAlgorithms< DataTypes >::Real vol=
            dot(tetraPos[1]-tetraPos[0],cross(tetraPos[2]-tetraPos[0],tetraPos[3]-tetraPos[0]));
        size_t i;
		for (i=0;i<4;++i) {
			if (i%2==0) 
				v[i]=dot(tetraPos[(i+1)%4]-pos,cross(tetraPos[(i+2)%4]-pos,tetraPos[(i+3)%4]-pos))/vol;
			else
				v[i]=-dot(tetraPos[(i+1)%4]-pos,cross(tetraPos[(i+2)%4]-pos,tetraPos[(i+3)%4]-pos))/vol;

		}
		// check it sums to 1
		assert(fabs(v[0]+v[1]+v[2]+v[3]-1)<1e-6);
		return(v);
}

template<class DataTypes>
void AdaptiveBezierTetrahedronSetTopologyAlgorithms< DataTypes >::updateTetrahedronDegree(
	sofa::helper::vector<TetraID>& loweringDegreeTetrahedra,
	sofa::helper::vector<TetraID>& raisingDegreeTetrahedra) 
{
	// points to be removed in descending order
	sofa::helper::set<size_t,std::greater<size_t> > pointsToBeRemoved;
	HighOrderDegreeType degree=m_adaptiveContainer->getDegree();
	size_t i,j,p;
	size_t tetrahedronIndex;
	sofa::helper::vector<unsigned int> bezierTetrahedronChangedArray;
	sofa::helper::vector<unsigned int> bezierTetrahedronChangedTriangleArray;
	sofa::helper::vector<unsigned int> bezierTetrahedronChangedEdgeArray;
	std::set<unsigned int> bezierTetrahedronChangedSet;

	helper::WriteOnlyAccessor<Data<SeqBezierDegree> >  tetrahedronDegree=m_adaptiveContainer->d_tetrahedronDegree;
	helper::WriteOnlyAccessor<Data<SeqBezierDegree> >  edgeDegree=m_adaptiveContainer->d_edgeDegree;
	helper::WriteOnlyAccessor<Data<SeqBezierDegree> >  triangleDegree=m_adaptiveContainer->d_triangleDegree;

	const VecCoord& coords =(this->object->read(core::ConstVecCoordId::position())->getValue());

	bool useSurfaceExtrapolation=m_adaptiveContainer->d_useSurfaceExtrapolation.getValue();

	size_t pointsPerTriangle=(degree-1)*(degree-2)/2;
	size_t pointsPerTetrahedron=(degree-1)*(degree-2)*(degree-3)/6;

	// process the tetrahedra that need to be transformed as linear
	if (!loweringDegreeTetrahedra.empty()) {
		Loc2GlobMapIterator itl2g;
		ControlPointLocation cpl;
		size_t k;
		Coord pos,tetraPos[4];
		
		for (p=0;p<loweringDegreeTetrahedra.size();++p) {
			assert(tetrahedronDegree[loweringDegreeTetrahedra[p]]==degree);
			tetrahedronIndex=loweringDegreeTetrahedra[p];
			 
			if (degree>1) {
				EdgesInTetrahedron eit=m_adaptiveContainer->getEdgesInTetrahedron(tetrahedronIndex);
				std::map<size_t, std::pair<size_t,size_t> >::iterator itedge; 


				for (i=0;i<6;++i) {
					Edge e=m_adaptiveContainer->getEdge(eit[i]);
					bool onSurface=false;
					// if the edge on the surface then store baycentric coordinates
					if ((itedge=m_adaptiveContainer->edgeOnSurfaceMap.find(eit[i]))!=m_adaptiveContainer->edgeOnSurfaceMap.end()) {
						onSurface=true;
					}
					// test if the edge is linear or of full degree
					if (edgeDegree[eit[i]]==degree) {
						bezierTetrahedronChangedEdgeArray.push_back(eit[i]);
						// add the tetrahedra adjacent to that edge in the set of tetrahedra to be updated
						const BaseMeshTopology::TetrahedraAroundEdge&tae= m_adaptiveContainer->getTetrahedraAroundEdge(eit[i]);
						for(j=0;j<tae.size();++j) {
							bezierTetrahedronChangedSet.insert(tae[j]);
						}
						Edge e=m_adaptiveContainer->getEdge(eit[i]);
						for(j=0;j<(degree-1);++j) {
							cpl=ControlPointLocation(eit[i],std::make_pair(HighOrderTetrahedronSetTopologyContainer::EDGE,j));
							// get the index of the associated DOF from the map.
							itl2g=m_adaptiveContainer->locationToGlobalIndexMap.find(cpl);
							// check that the DOF exists.
							assert(itl2g!=m_adaptiveContainer->locationToGlobalIndexMap.end());
							// add dof to the list
							pointsToBeRemoved.insert((*itl2g).second);
							if (onSurface) {
								// store barycentric coordinates
								// dof position
								pos=coords[(*itl2g).second];
								if (useSurfaceExtrapolation==false) {
									// get the 2 indices of adjacent tetrahedra
									std::map<size_t, std::pair<size_t,size_t> >::iterator ite1=m_adaptiveContainer->edgeOnSurfaceMap.find(eit[i]);
									assert(ite1!=m_adaptiveContainer->edgeOnSurfaceMap.end());
									Tetrahedron tet1=m_adaptiveContainer->getTetrahedron((*ite1).second.first);
									for(k=0;k<4;++k) {
										tetraPos[k]=coords[tet1[k]];
									}
									Vec4 bar1=getBarycentricCoordinates<DataTypes>(pos,tetraPos)/2.0f;
									Tetrahedron tet2=m_adaptiveContainer->getTetrahedron((*ite1).second.second);
									for(k=0;k<4;++k) {
										tetraPos[k]=coords[tet2[k]];
									}
									Vec4 bar2=getBarycentricCoordinates<DataTypes>(pos,tetraPos)/2.0f;
									std::pair<size_t,size_t> edgeLocation=std::make_pair(eit[i],j);
									// store the pair of barycentric coordinates
									m_adaptiveContainer->edgeControlPointsBarycentricCoord[edgeLocation]=std::make_pair(bar1,bar2);
								} else {
									// use surface displacement based extrapolation
									// get rest position of current point
									Coord resPosition=(*(restPositionMap.find(cpl))).second;
									// get rest position of the 2 extremities
									Real w=(Real)(degree-j-1)/(Real)degree;
									Coord offset=pos-resPosition-
										(coords[e[0]]-restPositionArray[e[0]])*w -
										(coords[e[1]]-restPositionArray[e[1]])*(1-w);
									offsetPositionMap[cpl]=offset;
								}
							}
						}

						// set degree to linear
						edgeDegree[eit[i]]=1;
					}
				}
			}
			// triangle index
			if (degree>2) {
				TrianglesInTetrahedron tit=m_adaptiveContainer->getTrianglesInTetrahedron(tetrahedronIndex);

				for (i=0;i<4;++i) {
					// test if the triangle is linear or of degree d
					if (triangleDegree[tit[i]]==degree) {
						bezierTetrahedronChangedTriangleArray.push_back(tit[i]);
						std::map<size_t,size_t>::iterator ittriangle; 
						bool onSurface=false;
						// if the edge on the surface then store baycentric coordinates
						if ((ittriangle=m_adaptiveContainer->trianglesOnSurfaceMap.find(tit[i]))!=m_adaptiveContainer->trianglesOnSurfaceMap.end()) {
							onSurface=true;
						}
						// add the 2 tetrahedra adjacent to that triangle in the set of tetrahedra to be updated
						const BaseMeshTopology::TetrahedraAroundTriangle&tat= m_adaptiveContainer->getTetrahedraAroundTriangle(tit[i]);
						for(j=0;j<tat.size();++j) {
							bezierTetrahedronChangedSet.insert(tat[j]);
						}

						for (j=0;j<pointsPerTriangle; ++j) {
							cpl=ControlPointLocation(tit[i],std::make_pair(HighOrderTetrahedronSetTopologyContainer::TRIANGLE,j));
							// get the index of the associated DOF from the map.
							itl2g=m_adaptiveContainer->locationToGlobalIndexMap.find(cpl);
							// check that the DOF exists.
							assert(itl2g!=m_adaptiveContainer->locationToGlobalIndexMap.end());
							// add dof to the list
							pointsToBeRemoved.insert((*itl2g).second);
							if (onSurface) {	
								Triangle tr=m_adaptiveContainer->getTriangle(tit[i]);
								// dof position
								pos=coords[(*itl2g).second];
								if (useSurfaceExtrapolation==false) {
									// store barycentric coordinates
									// get the index of adjacent tetrahedra
									std::map<size_t, size_t>::iterator itt1=m_adaptiveContainer->trianglesOnSurfaceMap.find(tit[i]);
									assert(itt1!=m_adaptiveContainer->trianglesOnSurfaceMap.end());
									Tetrahedron tet1=m_adaptiveContainer->getTetrahedron((*itt1).second);
									for(k=0;k<4;++k) {
										tetraPos[k]=coords[tet1[k]];
									}
									Vec4 bar=getBarycentricCoordinates<DataTypes>(pos,tetraPos);

									std::pair<size_t,size_t> triangleLocation=std::make_pair(tit[i],j);
									// store the pair of barycentric coordinates
									m_adaptiveContainer->triangleControlPointsBarycentricCoord[triangleLocation]=bar;
								} else {
									// use surface displacement based extrapolation
									// get rest position of current point
									Coord restPosition=(*(restPositionMap.find(cpl))).second;
									// get rest position of the 2 extremities
									TriangleIndexVector tbi;
									m_adaptiveContainer->getTriangleIndexFromTriangleOffset(cpl.second.second,tbi);
									Coord offset=pos-restPosition-
										(coords[tr[0]]-restPositionArray[tr[0]])*tbi[0]/(Real)degree-
										(coords[tr[1]]-restPositionArray[tr[1]])*tbi[1]/(Real)degree-
										(coords[tr[2]]-restPositionArray[tr[2]])*tbi[2]/(Real)degree;
									offsetPositionMap[cpl]=offset;
								}
							} 
						}
						
						// set degree to linear
						triangleDegree[tit[i]]=1;
					}
				}
			}
			// tetrahedron index
			if (degree>3) {

				bezierTetrahedronChangedSet.insert(tetrahedronIndex);

				for (j=0;j<pointsPerTetrahedron; ++j) {
					cpl=ControlPointLocation(tetrahedronIndex,std::make_pair(HighOrderTetrahedronSetTopologyContainer::TETRAHEDRON,j));
					// get the index of the associated DOF from the map.
					itl2g=m_adaptiveContainer->locationToGlobalIndexMap.find(cpl);
					// check that the DOF exists.
					assert(itl2g!=m_adaptiveContainer->locationToGlobalIndexMap.end());
					// add dof to the list
					pointsToBeRemoved.insert((*itl2g).second);
				}
			}

			// set tetrahedron degree to 1								
			tetrahedronDegree[tetrahedronIndex]=1;
		}


		if (!pointsToBeRemoved.empty()) {
			// now send the warning to topology containers.
			sofa::helper::vector<unsigned int> indexList(pointsToBeRemoved.size());
			std::copy(pointsToBeRemoved.begin(),pointsToBeRemoved.end(),indexList.begin());
			this->m_modifier->removePointsWarning(indexList);
			this->m_modifier->propagateTopologicalChanges();
			// update the weight array 
			sofa::helper::set<size_t>::iterator itset,itsetNext;
            unsigned int lastPoint = this->m_container->getNbPoints() - 1;
			HighOrderTetrahedronSetTopologyContainer::SeqWeights &wa=*(m_adaptiveContainer->d_weightArray.beginEdit());	
			for (itset=pointsToBeRemoved.begin();itset!=pointsToBeRemoved.end();++itset,--lastPoint)
			{
				if ((*itset)!=lastPoint)
					wa[(*itset)]=wa[lastPoint];
			}
			wa.resize(wa.size()-pointsToBeRemoved.size());
			m_adaptiveContainer->d_weightArray.endEdit();	
			// update the 2 maps.
			size_t offset;
			Glob2LocMapIterator itg2l,itg2l2,itg2lend; 
			// for each point to be removed in increasing order
			offset=1;
            lastPoint = this->m_container->getNbPoints() - 1;
			for (itset=pointsToBeRemoved.begin();itset!=pointsToBeRemoved.end();++itset,--lastPoint)
			{
				// swap the index of itset with the last entry point.

				itg2l=m_adaptiveContainer->globalIndexToLocationMap.find((*itset));
				assert(itg2l!=m_adaptiveContainer->globalIndexToLocationMap.end());
				itl2g=m_adaptiveContainer->locationToGlobalIndexMap.find((*itg2l).second);
				assert(itl2g!=m_adaptiveContainer->locationToGlobalIndexMap.end());
				m_adaptiveContainer->locationToGlobalIndexMap.erase(itl2g);
				m_adaptiveContainer->globalIndexToLocationMap.erase(itg2l);
				if (lastPoint!=(*itset)) {
					// swap the 2 indices
					itg2l=m_adaptiveContainer->globalIndexToLocationMap.find(lastPoint);
					assert(itg2l!=m_adaptiveContainer->globalIndexToLocationMap.end());
					ControlPointLocation cpl=(*itg2l).second;
					m_adaptiveContainer->globalIndexToLocationMap.erase(itg2l);
					itl2g=m_adaptiveContainer->locationToGlobalIndexMap.find(cpl);
					assert(itl2g!=m_adaptiveContainer->locationToGlobalIndexMap.end());
					m_adaptiveContainer->locationToGlobalIndexMap[cpl]=(*itset);
					m_adaptiveContainer->globalIndexToLocationMap[*itset]=cpl;
					// update the tetrahedronDOFArray
					// find the list of tetrahedra where the index should be swaped
					std::vector<size_t> tetrahedronList;
					
					if (cpl.second.first==HighOrderTetrahedronSetTopologyContainer::TETRAHEDRON) {
						tetrahedronList.push_back(cpl.first);
					} else if (cpl.second.first==HighOrderTetrahedronSetTopologyContainer::TRIANGLE) {
						const TetrahedraAroundTriangle &tat=m_adaptiveContainer->getTetrahedraAroundTriangle(cpl.first);
						for (j=0;j<tat.size();++j) {
							tetrahedronList.push_back(tat[j]);
						}
					}  else if (cpl.second.first==HighOrderTetrahedronSetTopologyContainer::EDGE) {
						const TetrahedraAroundEdge &tae=m_adaptiveContainer->getTetrahedraAroundEdge(cpl.first);
						for (j=0;j<tae.size();++j) {
							tetrahedronList.push_back(tae[j]);
						}
					} else {
						assert(false);
					}
					// swap indices
					AdaptiveHighOrderTetrahedronSetTopologyContainer::BezierDOFInTetrahedron da,daCopy;
					for (j=0;j<tetrahedronList.size();++j) {
						// find the DOF index associated with lastPoint
						da=m_adaptiveContainer->tetrahedronDOFArray[tetrahedronList[j]];
						for (k=0;da[k]!=lastPoint;++k);
						assert(k<da.size());
						da[k]=(*itset);
						m_adaptiveContainer->tetrahedronDOFArray[tetrahedronList[j]]=da;
					}
				}
			}
			
            lastPoint = this->m_container->getNbPoints() - 1;
			for (itset=pointsToBeRemoved.begin();itset!=pointsToBeRemoved.end();++itset,--lastPoint)
			{

			}
			// actually remove those points
			this->m_modifier->removePointsProcess(indexList);
			this->m_modifier->propagateTopologicalChanges();

		}


	}
	// process the tetrahedra that need to be degree raised from linear
	if (!raisingDegreeTetrahedra.empty()) {
		// first add points
		// collect information about control points to be added
		std::vector<Coord> pointsToBeAdded;
		// for each new control points defines its position as the barycentric coordinates from ancestors
		sofa::helper::vector< sofa::helper::vector<unsigned int> > ancestors;
		sofa::helper::vector< sofa::helper::vector<double> > coefs;

		HighOrderTetrahedronSetTopologyContainer::SeqWeights &wa=*(m_adaptiveContainer->d_weightArray.beginEdit());	

		Loc2GlobMapIterator itl2g;
		ControlPointLocation cpl;
		size_t k;


		// the global index of the first added point 
		size_t initIndex=m_adaptiveContainer->getNbPoints();
		size_t rank=initIndex;
		assert((m_adaptiveContainer->globalIndexToLocationMap.size())==initIndex);

		// parse all tetrahedron to be degree raised
		for (p=0;p<raisingDegreeTetrahedra.size();++p) {
			assert(tetrahedronDegree[raisingDegreeTetrahedra[p]]==1);
			tetrahedronIndex=raisingDegreeTetrahedra[p];
			Tetrahedron tet=m_adaptiveContainer->getTetrahedron(tetrahedronIndex);
			// add edge points

			EdgesInTetrahedron eit=m_adaptiveContainer->getEdgesInTetrahedron(tetrahedronIndex);
			for (i=0;i<6;++i) {
				Edge e=m_adaptiveContainer->getEdge(eit[i]);
				// test if the edge is linear or of full degree
				if(edgeDegree[eit[i]]==1) {
					bezierTetrahedronChangedEdgeArray.push_back(eit[i]);
					// add the tetrahedra adjacent to that edge in the set of tetrahedra to be updated
					const BaseMeshTopology::TetrahedraAroundEdge&tae= m_adaptiveContainer->getTetrahedraAroundEdge(eit[i]);
					for(j=0;j<tae.size();++j) {
						bezierTetrahedronChangedSet.insert(tae[j]);
					}

					std::map<size_t, std::pair<size_t,size_t> >::iterator itedge; 
					bool onSurface=false;
					
					sofa::helper::vector<unsigned int>  ancestorArray;
					
					// if the edge on the surface then store baycentric coordinates
					if ((itedge=m_adaptiveContainer->edgeOnSurfaceMap.find(eit[i]))!=m_adaptiveContainer->edgeOnSurfaceMap.end()) {
						onSurface=true;
						// get the 2 adjacent tetrahedra of the edge
						std::map<size_t, std::pair<size_t,size_t> >::iterator  itet;
						itet=m_adaptiveContainer->edgeOnSurfaceMap.find(eit[i]);
						assert(itet!=m_adaptiveContainer->edgeOnSurfaceMap.end());
						Tetrahedron tet1=m_adaptiveContainer->getTetrahedron((*itet).second.first);
						Tetrahedron tet2=m_adaptiveContainer->getTetrahedron((*itet).second.second);
						// get the indices of the 2 tetrahedra
						for (k=0;k<4;++k)
							ancestorArray.push_back(tet1[k]);
						for (k=0;k<4;++k)
							ancestorArray.push_back(tet2[k]);

					} else {
						for (k=0;k<2;++k) 
							ancestorArray.push_back(e[k]);
					}
					for(j=0;j<(degree-1);++j) {
						// define control points as barycentric coordinates from edge extremities
						
						sofa::helper::vector<double>  coefArray;

						cpl=ControlPointLocation(eit[i],std::make_pair(HighOrderTetrahedronSetTopologyContainer::EDGE,j));
						wa.push_back(m_adaptiveContainer->bezierWeightsMap[cpl]);

						if (onSurface) {
							if (useSurfaceExtrapolation==false) {
								std::pair<size_t,size_t> edgeLocation=std::make_pair(eit[i],j);
								std::map<std::pair<size_t,size_t>,std::pair<Vec4,Vec4> >::iterator itew;
								// get the previously barycentric coordinates of the edge points
								itew=m_adaptiveContainer->edgeControlPointsBarycentricCoord.find(edgeLocation);
								assert(itew!=m_adaptiveContainer->edgeControlPointsBarycentricCoord.end());

								for(k=0;k<4;++k) {
									coefArray.push_back((*itew).second.first[k]);
								}
								for(k=0;k<4;++k) {
									coefArray.push_back((*itew).second.second[k]);
								}
							} else {
								// in this case the position of the vertex is computed from its rest position and the interpolated
								// displacement field

								Coord restPosition=restPositionMap[cpl];
								Edge e=m_adaptiveContainer->getEdge(eit[i]);
								Coord pos=restPosition+//offsetPositionMap[cpl]+
									(coords[e[0]]-restPositionArray[e[0]])*(degree-j-1)/(Real)degree +
									(coords[e[1]]-restPositionArray[e[1]])*(j+1)/(Real)degree;
								// now must provide the barycentric coordinate of that point
								// with respect to the 2 adjacent tetrahedra
								Coord tetraPos[4];
								for (k=0;k<4;++k) 
									tetraPos[k]=coords[ancestorArray[k]];
								Vec4 bar1=getBarycentricCoordinates<DataTypes>(pos,tetraPos)/2.0f;
								for (k=0;k<4;++k) 
									tetraPos[k]=coords[ancestorArray[k+4]];
								Vec4 bar2=getBarycentricCoordinates<DataTypes>(pos,tetraPos)/2.0f;
								for(k=0;k<4;++k) 
									coefArray.push_back(bar1[k]);
								for(k=0;k<4;++k) 
									coefArray.push_back(bar2[k]);
							}
						} else {
							coefArray.push_back((Real)(degree-j-1)/(Real)degree);
							coefArray.push_back((Real)(j+1)/(Real)degree);
						}
						coefs.push_back(coefArray);
						ancestors.push_back(ancestorArray);
						// update the 2 maps
						m_adaptiveContainer->locationToGlobalIndexMap.insert(
							std::pair<ControlPointLocation,size_t>(cpl,rank));
						m_adaptiveContainer->globalIndexToLocationMap[rank]=cpl;
						rank++;
					}
					// set degree to full degree
					edgeDegree[eit[i]]=degree;
				}
			}

			// triangle index
			size_t l,m;
			TrianglesInTetrahedron tit=m_adaptiveContainer->getTrianglesInTetrahedron(tetrahedronIndex);
			for (i=0;i<4;++i) {
				// verify that the triangle is linear
				if (triangleDegree[tit[i]]==1) {
					bezierTetrahedronChangedTriangleArray.push_back(tit[i]);
					// add the 2 tetrahedra adjacent to that triangle in the set of tetrahedra to be updated
					const BaseMeshTopology::TetrahedraAroundTriangle&tat= m_adaptiveContainer->getTetrahedraAroundTriangle(tit[i]);
					for(j=0;j<tat.size();++j) {
						bezierTetrahedronChangedSet.insert(tat[j]);
					}

					sofa::helper::vector<unsigned int>  ancestorArray;
					std::map<size_t,size_t>::iterator ittriangle; 
					bool onSurface=false;
					// if the edge on the surface then store baycentric coordinates
					if ((ittriangle=m_adaptiveContainer->trianglesOnSurfaceMap.find(tit[i]))!=m_adaptiveContainer->trianglesOnSurfaceMap.end()) {
						onSurface=true;

						// get the  adjacent tetrahedron of the triangle
						std::map<size_t, size_t >::iterator  itet;
						itet=m_adaptiveContainer->trianglesOnSurfaceMap.find(tit[i]);
						assert(itet!=m_adaptiveContainer->trianglesOnSurfaceMap.end());
						Tetrahedron tet1=m_adaptiveContainer->getTetrahedron((*itet).second);
						// get the indices of the adjacent tetrahedra
						for (k=0;k<4;++k)
							ancestorArray.push_back(tet1[k]);

					} else {
						Triangle tr=m_adaptiveContainer->getTriangle(tit[i]);
						for (l=0;l<3;++l) 
							ancestorArray.push_back(tr[l]);
					}
					for (m=0,j=1;j<(size_t)(degree-1);++j) {
						for (k=1;k<(degree-j);++k,++m) {
							// update weights
							cpl=ControlPointLocation(tit[i],std::make_pair(HighOrderTetrahedronSetTopologyContainer::TRIANGLE,m));
							wa.push_back(m_adaptiveContainer->bezierWeightsMap[cpl]);

							// define control points as barycentric coordinates from edge extremities

							sofa::helper::vector<double>  coefArray;
							if (onSurface) {
								if (useSurfaceExtrapolation==false) {
									std::pair<size_t,size_t> triangleLocation=std::make_pair(tit[i],m);
									std::map<std::pair<size_t,size_t>,Vec4 >::iterator itew;
									// get the previously barycentric coordinates of the triangle points
									itew=m_adaptiveContainer->triangleControlPointsBarycentricCoord.find(triangleLocation);
									assert(itew!=m_adaptiveContainer->triangleControlPointsBarycentricCoord.end());

									for(l=0;l<4;++l) {
										coefArray.push_back((*itew).second[l]);
									}
								} else {
									// in this case the position of the vertex is computed from its rest position and the interpolated
									// displacement field
																	Coord restPosition=restPositionMap[cpl];
									Triangle  tr=m_adaptiveContainer->getTriangle(tit[i]);
									Coord pos=restPosition+offsetPositionMap[cpl]+
										(coords[tr[0]]-restPositionArray[tr[0]])*(j)/(Real)degree +
										(coords[tr[1]]-restPositionArray[tr[1]])*(k)/(Real)degree +
										(coords[tr[2]]-restPositionArray[tr[2]])*(degree-j-k)/(Real)degree;
									// now must provide the barycentric coordinate of that point
									// with respect to the 2 adjacent tetrahedra
									Coord tetraPos[4];
									for (l=0;l<4;++l) 
										tetraPos[l]=coords[ancestorArray[l]];
									Vec4 bar1=getBarycentricCoordinates<DataTypes>(pos,tetraPos);
									for (l=0;l<4;++l) 
										coefArray.push_back(bar1[l]);
								}
							} else {
								coefArray.push_back((Real)(j)/(Real)degree);
								coefArray.push_back((Real)(k)/(Real)degree);
								coefArray.push_back((Real)(degree-j-k)/(Real)degree);
							}


							coefs.push_back(coefArray);
							ancestors.push_back(ancestorArray);
							// update the 2 maps
							m_adaptiveContainer->locationToGlobalIndexMap.insert(
								std::pair<ControlPointLocation,size_t>(cpl,rank));
							m_adaptiveContainer->globalIndexToLocationMap[rank]=cpl;
							rank++;

						}
					}
					// set degree to full degree
					triangleDegree[tit[i]]=degree;

				}
			}
			// tetrahedron index

			// verify that the tetrahedra is linear
			assert(tetrahedronDegree[tetrahedronIndex]==1);
			// test if the triangle is linear or of degree d
			if (degree>3) {
				bezierTetrahedronChangedSet.insert(tetrahedronIndex);
				sofa::helper::vector<unsigned int>  ancestorArray;
				for (l=0;l<4;++l) 
					ancestorArray.push_back(tet[l]);
				for (m=0,i=1;i<(size_t)(degree-2);++i) {
					for (j=1;j<(degree-i-1);++j) {
						for (k=1;k<(degree-j-i);++k,++m) {

							cpl=ControlPointLocation(tetrahedronIndex,std::make_pair(HighOrderTetrahedronSetTopologyContainer::TETRAHEDRON,m));
							wa.push_back(m_adaptiveContainer->bezierWeightsMap[cpl]);

							// define control points as barycentric coordinates from edge extremities
							
							sofa::helper::vector<double>  coefArray;
							coefArray.push_back((Real)(i)/(Real)degree);
							coefArray.push_back((Real)(j)/(Real)degree);
							coefArray.push_back((Real)(k)/(Real)degree);
							coefArray.push_back((Real)(degree-i-j-k)/(Real)degree);
							
							coefs.push_back(coefArray);
							ancestors.push_back(ancestorArray);
							// update the 2 maps
							m_adaptiveContainer->locationToGlobalIndexMap.insert(
								std::pair<ControlPointLocation,size_t>(cpl,rank));
							m_adaptiveContainer->globalIndexToLocationMap[rank]=cpl;
							rank++;
						}
					}
				}
			}
			// set degree to full degree
			tetrahedronDegree[tetrahedronIndex]=degree;
		}
		// actually add the points
        this->m_modifier->addPointsProcess(ancestors.size());
        this->m_modifier->addPointsWarning(ancestors.size(), ancestors, coefs, true);
		this->m_modifier->propagateTopologicalChanges();
		}
		// now parse the set of tetrahedra to be updated
		std::set<unsigned int>::iterator itet=bezierTetrahedronChangedSet.begin();
		AdaptiveHighOrderTetrahedronSetTopologyContainer::SeqWeightedDOFArray swda;
		AdaptiveHighOrderTetrahedronSetTopologyContainer::BezierDOFInTetrahedron da;
		for (;itet!=bezierTetrahedronChangedSet.end();++itet) {
			bezierTetrahedronChangedArray.push_back(*itet);
			// now update entry in m_adaptiveContainer about each tetrahedron
			swda.clear();
			da.clear();
			m_adaptiveContainer->getWeightedDOFOfTetrahedron((*itet),swda,da);
			m_adaptiveContainer->weightedDOFArray[(*itet)]=swda;
			this->m_adaptiveContainer->tetrahedronDOFArray[(*itet)]=da;
		}
#ifndef NDEBUG
		// check consistency between both maps
		assert(m_adaptiveContainer->locationToGlobalIndexMap.size()==m_adaptiveContainer->globalIndexToLocationMap.size());
		std::map<ControlPointLocation,size_t>::iterator itcpl;
		std::map<size_t,ControlPointLocation>::iterator itgi;
		for (itcpl=m_adaptiveContainer->locationToGlobalIndexMap.begin();
			itcpl!=m_adaptiveContainer->locationToGlobalIndexMap.end();++itcpl) {
			itgi=m_adaptiveContainer->globalIndexToLocationMap.find(itcpl->second);
			assert(itgi!=m_adaptiveContainer->globalIndexToLocationMap.end());
			assert((*itgi).second==(*itcpl).first);
		}
#endif
		// now throw topological event
        this->m_adaptiveModifier->changeDegree(bezierTetrahedronChangedEdgeArray,
			bezierTetrahedronChangedTriangleArray,
			bezierTetrahedronChangedArray);
        this->m_modifier->notifyEndingEvent();
        this->m_modifier->propagateTopologicalChanges();
#ifndef NDEBUG
		m_adaptiveContainer->checkTopology();
#endif

}


} // namespace topology

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENTS_TETEAHEDRONSETTOPOLOGYALGORITHMS_INL
