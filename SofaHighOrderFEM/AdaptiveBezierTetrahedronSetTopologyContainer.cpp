
#include "AdaptiveBezierTetrahedronSetTopologyContainer.h"
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>
//#include <tuple>
//#include "boost/tuple/tuple.hpp"

namespace sofa
{

namespace component
{

namespace topology
{

using namespace std;
using namespace sofa::defaulttype;

SOFA_DECL_CLASS(AdaptiveHighOrderTetrahedronSetTopologyContainer)
int AdaptiveHighOrderTetrahedronSetTopologyContainerClass = core::RegisterObject("Adaptive High order Tetrahedron set topology container")
        .add< AdaptiveHighOrderTetrahedronSetTopologyContainer >()
        ;

const unsigned int edgesInTetrahedronArray[6][2] = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}};
///convention triangles in tetra (orientation interior)
const unsigned int trianglesInTetrahedronArray[4][3]= {{1,2,3}, {0,3,2}, {1,3,0},{0,2,1}};

typedef std::map<AdaptiveHighOrderTetrahedronSetTopologyContainer::ControlPointLocation,size_t>::iterator Loc2GlobMapIterator;
typedef std::map<TetrahedronIndexVector,size_t>::iterator OffsetMapIterator;

AdaptiveHighOrderTetrahedronSetTopologyContainer::AdaptiveHighOrderTetrahedronSetTopologyContainer()
    : HighOrderTetrahedronSetTopologyContainer()
	, d_tetrahedronDegree(initData(&d_tetrahedronDegree, SeqBezierDegree(),"tetrahedronDegree", "Gives the degree of each Bezier tetrahedron "))
	, d_triangleDegree(initData(&d_triangleDegree, SeqBezierDegree(),"triangleDegree", "Gives the degree of each Bezier triangle "))
	, d_edgeDegree(initData(&d_edgeDegree, SeqBezierDegree(),"edgeDegree", "Gives the degree of each Bezier edge "))
	, d_useSurfaceExtrapolation(initData(&d_useSurfaceExtrapolation, true,"useSurfaceExtrapolation", "indicates if surface or volumetric extrapolation is use to update the surface triangles"))

{
    
}



void AdaptiveHighOrderTetrahedronSetTopologyContainer::init()
{
	 d_tetrahedronDegree.updateIfDirty();
 	 d_triangleDegree.updateIfDirty();
  	 d_edgeDegree.updateIfDirty();

	HighOrderTetrahedronSetTopologyContainer::init(); // initialize the tetrahedron array
}
// returns a pointer to an array describing the DOFs of all tetrahedra
const AdaptiveHighOrderTetrahedronSetTopologyContainer::SeqSeqWeightedDOFArray &AdaptiveHighOrderTetrahedronSetTopologyContainer::getWeightedDOFArrayArray() const
{
	return weightedDOFArray;
}
// returns a pointer to an array describing the DOFs of all tetrahedra
const AdaptiveHighOrderTetrahedronSetTopologyContainer::SeqBezierDOFInTetrahedron &AdaptiveHighOrderTetrahedronSetTopologyContainer::getTetrahedronDOFArrayArray() const {
	return tetrahedronDOFArray;
}
// returns a pointer to an array describing the DOFs of all tetrahedra
const AdaptiveHighOrderTetrahedronSetTopologyContainer::SeqWeightedDOFArray 
	&AdaptiveHighOrderTetrahedronSetTopologyContainer::getWeightedDOFArray(const unsigned int i) {
		return weightedDOFArray[i];
}
// returns a pointer to an array describing the DOFs of all tetrahedra
const AdaptiveHighOrderTetrahedronSetTopologyContainer::BezierDOFInTetrahedron 
	&AdaptiveHighOrderTetrahedronSetTopologyContainer::getTetrahedronDOFArray(const unsigned int i){
		return tetrahedronDOFArray[i];
}
void AdaptiveHighOrderTetrahedronSetTopologyContainer::reinit()
{
	HighOrderTetrahedronSetTopologyContainer::reinit();
	// fill the 2 maps
	size_t i,j,rank;
	for (rank=0;rank<getNumberOfTetrahedralPoints();++rank) {
		ControlPointLocation cpl(rank,std::make_pair(HighOrderTetrahedronSetTopologyContainer::POINT,0));
		locationToGlobalIndexMap.insert(std::pair<ControlPointLocation,size_t>(cpl,rank));
		globalIndexToLocationMap[rank]=cpl;
	}
	size_t degree=getDegree();

	if (degree>1) {
		for (i=0;i<getNbEdges();++i) {
			for (j=1;j<degree;++j,rank++) {
				ControlPointLocation cpl(i,std::make_pair(HighOrderTetrahedronSetTopologyContainer::EDGE,j-1));
				locationToGlobalIndexMap.insert(std::pair<ControlPointLocation,size_t>(cpl,rank));
				globalIndexToLocationMap[rank]=cpl;
			}
		}
	}
	if (degree>2) {
		size_t k,l;
		
		for (i=0;i<getNbTriangles();++i) {
			for (l=0,j=1;j<(size_t)(degree-1);++j) {
				for (k=1;k<(degree-j);++k,++rank,++l) {
					ControlPointLocation cpl(i,std::make_pair(HighOrderTetrahedronSetTopologyContainer::TRIANGLE,l));
					locationToGlobalIndexMap.insert(std::pair<ControlPointLocation,size_t>(cpl,rank));
					globalIndexToLocationMap[rank]=cpl;
				}
			}
		}
	}
	if (degree>3) {
		size_t k,l,m;
		
		for (i=0;i<getNbTetrahedra();++i) {
			for (m=0,j=1;j<(size_t)(degree-2);++j) {
				for (k=1;k<(size_t)(degree-j-1);++k) {
					for (l=1;l<(degree-j-k);++l,++rank,++m) {
						ControlPointLocation cpl(i,std::make_pair(HighOrderTetrahedronSetTopologyContainer::TETRAHEDRON,m));
						locationToGlobalIndexMap.insert(std::pair<ControlPointLocation,size_t>(cpl,rank));
						globalIndexToLocationMap[rank]=cpl;
					}
				}
			}
		}
	}
	// initialize degree arrays with full degree
	helper::WriteOnlyAccessor<Data<SeqBezierDegree> >  tetrahedronDegree=d_tetrahedronDegree;
	tetrahedronDegree.resize(getNbTetrahedra());
	std::fill(tetrahedronDegree.begin(),tetrahedronDegree.end(),degree);

	helper::WriteOnlyAccessor<Data<SeqBezierDegree> >  triangleDegree=d_triangleDegree;
	triangleDegree.resize(getNbTriangles());
	std::fill(triangleDegree.begin(),triangleDegree.end(),degree);
	
	helper::WriteOnlyAccessor<Data<SeqBezierDegree> >  edgeDegree=d_edgeDegree;
	edgeDegree.resize(getNbEdges());
	std::fill(edgeDegree.begin(),edgeDegree.end(),degree);

	// initialize 	weightedDOFArray
	SeqWeightedDOFArray swda;
	BezierDOFInTetrahedron da;
	for (i=0;i<getNbTetrahedra();++i) {
		getWeightedDOFOfTetrahedron(i,swda,da);
		weightedDOFArray.push_back(swda);
	}
	
	
	// now store the triangles and edges of the surface in 2 sets
	for (i=0;i<getNbTriangles();++i) {
		// detect triangles on the surface
		if (getTetrahedraAroundTriangle(i).size()==1) {
			// add index of the adjacent tetrahedra in the map
			triangleOnSurfaceSet.insert(i);
			EdgesInTriangle eit=getEdgesInTriangle(i);
			// store the associated tetrahedron for each edge
			for(j=0;j<3;++j) {
				edgeOnSurfaceSet.insert(eit[j]);
			}		
		}
	}
	// for each edge on surface store the indices of the 2 adjacent tetrahedra
	std::multimap<size_t,size_t>  edgeMap;
	// now store the triangles and edges of the surface
	for (i=0;i<getNbTriangles();++i) {
		// detect triangles on the surface
		if (getTetrahedraAroundTriangle(i).size()==1) {
			// add index of the adjacent tetrahedra in the map
			trianglesOnSurfaceMap[i]=getTetrahedraAroundTriangle(i)[0];

			EdgesInTriangle eit=getEdgesInTriangle(i);
			// store the associated tetrahedron for each edge
			for(j=0;j<3;++j) {
				edgeMap.insert(std::pair<size_t,size_t>((size_t)(eit[j]),
					getTetrahedraAroundTriangle(i)[0]));

			}		
		}
	}

	std::multimap<size_t,size_t>::iterator it1,it2;
	std::set<size_t>::iterator itset;
	size_t tet1,tet2;
	std::pair< std::multimap<size_t,size_t>::iterator , std::multimap<size_t,size_t>::iterator  > ite;
	for (itset=edgeOnSurfaceSet.begin();itset!=edgeOnSurfaceSet.end();++itset) {
		ite=edgeMap.equal_range((*itset));
		assert(ite.first!=ite.second);
		it1=ite.first;
		tet1=(*it1).second;
		it1++;
		assert(it1!=ite.second);
		tet2=(*it1).second;
		it1++;
		// each edge should be associated with 2 tetrahedra and only 2
		assert(it1==ite.second);
		// update edgeMap
		edgeOnSurfaceMap[(*itset)]=std::pair<size_t,size_t>(tet1,tet2);
	}
	// store the bezier weights in a map
	SReal weight;
	ControlPointLocation cpl;
	size_t elementIndex,elementOffset;
	BezierTetrahedronPointLocation location;
	for (i=getNumberOfTetrahedralPoints();i<getNbPoints();++i) {
		weight=getWeight(i);
		this->getLocationFromGlobalIndex(i,location,elementIndex,elementOffset);
		cpl=ControlPointLocation(elementIndex,std::make_pair(location,elementOffset));
		bezierWeightsMap[cpl]=weight;
	}

}
size_t AdaptiveHighOrderTetrahedronSetTopologyContainer::getTetrahedronDegree(const size_t i) const {
	return(d_tetrahedronDegree.getValue()[i]);
}
const  AdaptiveHighOrderTetrahedronSetTopologyContainer::SeqBezierDegree &AdaptiveHighOrderTetrahedronSetTopologyContainer::getTetrahedronDegreeArray() const
{
	return(d_tetrahedronDegree.getValue());
}
size_t AdaptiveHighOrderTetrahedronSetTopologyContainer::getTriangleDegree(const size_t i) const {
	return(d_triangleDegree.getValue()[i]);
}
const  AdaptiveHighOrderTetrahedronSetTopologyContainer::SeqBezierDegree &AdaptiveHighOrderTetrahedronSetTopologyContainer::getTriangleDegreeArray() const
{
	return(d_triangleDegree.getValue());
}
size_t AdaptiveHighOrderTetrahedronSetTopologyContainer::getEdgeDegree(const size_t i) const {
	return(d_edgeDegree.getValue()[i]);
}
const  AdaptiveHighOrderTetrahedronSetTopologyContainer::SeqBezierDegree &AdaptiveHighOrderTetrahedronSetTopologyContainer::getEdgeDegreeArray() const
{
	return(d_edgeDegree.getValue());
}
bool AdaptiveHighOrderTetrahedronSetTopologyContainer::checkHighOrderTetrahedronTopology() 
{
	if (this->HighOrderTetrahedronSetTopologyContainer::checkHighOrderTetrahedronTopology()) {
		size_t i,j;
		size_t nbDofs=getNbPoints();
		assert(nbDofs==d_weightArray.getValue().size());
		for (i=0;i<getNumberOfTetrahedra();++i) {
			const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=
				getTetrahedronDOFArray(i);
			const AdaptiveHighOrderTetrahedronSetTopologyContainer::SeqWeightedDOFArray &swda=
				getWeightedDOFArray(i);
			// check the indices in  indexArray
			for(j=0;j<indexArray.size();++j) {
				assert(indexArray[j]<nbDofs);

			}
		}
		return true;
	}
}
void AdaptiveHighOrderTetrahedronSetTopologyContainer::getWeightedDOFOfTetrahedron(const size_t tetrahedronIndex, 
																				AdaptiveHighOrderTetrahedronSetTopologyContainer::SeqWeightedDOFArray &swda,
																				AdaptiveHighOrderTetrahedronSetTopologyContainer::BezierDOFInTetrahedron &da) 
{
	WeightedDOFArray wda;
	da.clear();
	swda.clear();

	Tetrahedron tet=getTetrahedron(tetrahedronIndex);
	// vertex index
	size_t i,j,k;
	for (i=0;i<4;++i) {
		wda.clear();
		da.push_back(tet[i]);
		wda.push_back(WeightedDOF(i,(SReal)1.0f));
		swda.push_back(wda);
	}
	// edge index
	HighOrderDegreeType degree=d_degree.getValue();
	if (degree>1) {
		Loc2GlobMapIterator itl2g;
		const AdaptiveHighOrderTetrahedronSetTopologyContainer::SeqBezierDegree &edgeDegree=d_edgeDegree.getValue();
		EdgesInTetrahedron eit=getEdgesInTetrahedron(tetrahedronIndex);
		for (i=0;i<6;++i) {
			Edge e=getEdge(eit[i]);
			// test if the edge is linear or of degree d
			if (edgeDegree[eit[i]]==degree) {
				// degree d
				for (j=0;j<(size_t)(degree-1);++j) {
					ControlPointLocation cpl;
					// test orientation of edge in the tetrahedron
					if (e[0]==tet[edgesInTetrahedronArray[i][0]]) {
						cpl=ControlPointLocation(eit[i],std::make_pair(HighOrderTetrahedronSetTopologyContainer::EDGE,j));
					} else {
						cpl=ControlPointLocation(eit[i],std::make_pair(HighOrderTetrahedronSetTopologyContainer::EDGE,degree-2-j));
					}
					// get the index of the associated DOF from the map.
					itl2g=locationToGlobalIndexMap.find(cpl);
					// check that the DOF exists.
					assert(itl2g!=locationToGlobalIndexMap.end());
					// set the mapped dof
					wda.clear();
					
					wda.push_back(WeightedDOF(da.size(),(SReal)1.0f));
					da.push_back((*itl2g).second);
					swda.push_back(wda);
				}

			} else {
				// linear  edge
				for (j=1;j<(size_t)(degree);++j) {
					// the DOF are linear mapped from the 2 extremities of the edge
					wda.clear();
					for(k=0;e[0]!=tet[k];++k);
						wda.push_back(WeightedDOF(k,(SReal)(degree-j)/degree));
		//			wda.push_back(WeightedDOF(k,(SReal)j/degree));
					for(k=0;e[1]!=tet[k];++k);
					wda.push_back(WeightedDOF(k,(SReal)j/degree));
	//				wda.push_back(WeightedDOF(k,(SReal)(degree-j)/degree));
					swda.push_back(wda);
				}
			}
			
		}
	}
	// triangle index
	if (degree>2) {
		TrianglesInTetrahedron tit=getTrianglesInTetrahedron(tetrahedronIndex);
		Loc2GlobMapIterator itl2g;
		const AdaptiveHighOrderTetrahedronSetTopologyContainer::SeqBezierDegree &triangleDegree=d_triangleDegree.getValue();
		size_t l;
		for (i=0;i<4;++i) {
			Triangle tr=getTriangle(tit[i]);
			// renumber the triangle as to fit the way it is stored in the reference Bezier tetrahedron
			Triangle indexTriangle;
			for (k=0;(tr[0]!=tet[trianglesInTetrahedronArray[i][k]]);++k);
			indexTriangle[k]=0;
			if (tr[1]==tet[trianglesInTetrahedronArray[i][(k+1)%3]]) {
				indexTriangle[(k+1)%3]=1;
				indexTriangle[(k+2)%3]=2;
			} else {
				indexTriangle[(k+2)%3]=1;
				indexTriangle[(k+1)%3]=2;
			}

			// test if the triangle is linear or of degree d
			if (triangleDegree[tit[i]]==degree) {
				ControlPointLocation cpl;
				for (j=1;j<(size_t)(degree-1);++j) {
					for (k=1;k<(degree-j);++k) {
						TetrahedronIndexVector bti(0,0,0,0);
						bti[trianglesInTetrahedronArray[i][indexTriangle[0]]]=j;
						bti[trianglesInTetrahedronArray[i][indexTriangle[1]]]=k;
						bti[trianglesInTetrahedronArray[i][indexTriangle[2]]]=degree-j-k;
						OffsetMapIterator omi=triangleOffsetMap.find(bti);
						cpl=ControlPointLocation(tit[i],std::make_pair(HighOrderTetrahedronSetTopologyContainer::TRIANGLE,(*omi).second));
						// get the index of the associated DOF from the map.
						itl2g=locationToGlobalIndexMap.find(cpl);
						// check that the DOF exists.
						assert(itl2g!=locationToGlobalIndexMap.end());
						// set the mapped dof
						wda.clear();
						//						wda.push_back(WeightedDOF((*itl2g).second,(SReal)1.0f));
						wda.push_back(WeightedDOF(da.size(),(SReal)1.0f));
						da.push_back((*itl2g).second);
						swda.push_back(wda);
					}
				}
			} else {
				for (j=1;j<(size_t)(degree-1);++j) {
					for (k=1;k<(degree-j);++k) {
						// the DOF are linear mapped from the 3 extremities of the triangle
						wda.clear();
				//		for(l=0;tr[0]!=tet[l];++l);
						wda.push_back(WeightedDOF(trianglesInTetrahedronArray[i][0],(SReal)j/degree));
				//		for(l=0;tr[1]!=tet[l];++l);
						wda.push_back(WeightedDOF(trianglesInTetrahedronArray[i][1],(SReal)k/degree));
						for(l=0;tr[2]!=tet[l];++l);
						wda.push_back(WeightedDOF(trianglesInTetrahedronArray[i][2],(SReal)(degree-j-k)/degree));
						swda.push_back(wda);
					}
				}
			}
		}
	}
	

	// tetrahedron index
	if (degree>3) {
		const AdaptiveHighOrderTetrahedronSetTopologyContainer::SeqBezierDegree &tetrahedronDegree=
			d_tetrahedronDegree.getValue();
		Loc2GlobMapIterator itl2g;
		// test if the triangle is linear or of degree d
		if (tetrahedronDegree[tetrahedronIndex]==degree) {
			ControlPointLocation cpl;
			size_t rank=0;
			for (i=1;i<(size_t)(degree-2);++i) {
				for (j=1;j<(degree-i-1);++j) {
					for (k=1;k<(degree-j-i);++k,++rank) {
						cpl=ControlPointLocation(tetrahedronIndex,std::make_pair(HighOrderTetrahedronSetTopologyContainer::TETRAHEDRON,rank));
						// get the index of the associated DOF from the map.
						itl2g=locationToGlobalIndexMap.find(cpl);
						// check that the DOF exists.
						assert(itl2g!=locationToGlobalIndexMap.end());
						// set the mapped dof
						wda.clear();
//						wda.push_back(WeightedDOF((*itl2g).second,(SReal)1.0f));
						wda.push_back(WeightedDOF(da.size(),(SReal)1.0f));
						da.push_back((*itl2g).second);
						swda.push_back(wda);
					}
				}
			}
		} else {
			for (i=1;i<(size_t)(degree-2);++i) {
				for (j=1;j<(degree-i-1);++j) {
					for (k=1;k<(degree-j-i);++k) {
						// the DOF are linear mapped from the 3 extremities of the triangle
						wda.clear();
						wda.push_back(WeightedDOF(0,(SReal)i/degree));
						wda.push_back(WeightedDOF(1,(SReal)j/degree));
						wda.push_back(WeightedDOF(2,(SReal)k/degree));
						wda.push_back(WeightedDOF(3,(SReal)(degree-j-k-i)/degree));
						swda.push_back(wda);
					}
				}
			}
		}
		
		
	}
	
}


} // namespace topology

} // namespace component

} // namespace sofa
