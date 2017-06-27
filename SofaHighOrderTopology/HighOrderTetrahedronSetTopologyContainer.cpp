

#include "HighOrderTetrahedronSetTopologyContainer.h"
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>
 

namespace sofa
{

namespace component
{

namespace topology
{

using namespace std;
using namespace sofa::defaulttype;

SOFA_DECL_CLASS(HighOrderTetrahedronSetTopologyContainer)
int HighOrderTetrahedronSetTopologyContainerClass = core::RegisterObject("High Order Tetrahedron set topology container")
        .add< HighOrderTetrahedronSetTopologyContainer >()
        ;

const unsigned int edgesInTetrahedronArray[6][2] = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}};
///convention triangles in tetra (orientation interior)
const unsigned int trianglesInTetrahedronArray[4][3]= {{1,2,3}, {0,3,2}, {1,3,0},{0,2,1}};



typedef std::pair<TetrahedronIndexVector,HighOrderTetrahedronSetTopologyContainer::ElementTetrahedronIndex> ElementMapType;
typedef std::map<TetrahedronIndexVector,HighOrderTetrahedronSetTopologyContainer::ElementTetrahedronIndex>::iterator ElementMapIterator;

typedef std::pair<TetrahedronIndexVector,size_t> OffsetMapType;
typedef std::map<TetrahedronIndexVector,size_t>::iterator OffsetMapIterator;
typedef std::map<TetrahedronIndexVector,size_t>::const_iterator OffsetMapConstIterator;

HighOrderTetrahedronSetTopologyContainer::HighOrderTetrahedronSetTopologyContainer()
    : TetrahedronSetTopologyContainer()
    , d_degree(initData(&d_degree, (size_t)0,"degree", "Degree of High Order Tetrahedra"))
    , d_numberOfTetrahedralPoints(initData(&d_numberOfTetrahedralPoints, (size_t) 0,"NbTetrahedralVertices", "Number of Tetrahedral Vertices"))
	, d_isRationalSpline(initData(&d_isRationalSpline, SeqBools(),"isRational", "If a bezier tetrahedron is rational or integral"))
	, d_weightArray(initData(&d_weightArray, SeqWeights(),"weights", "Array of weights for rational bezier tetrahedra"))
    , inputEdges(initData(&inputEdges,"inputHighOrderEdges","Edges where high order points are lying"))
    , inputTriangles(initData(&inputTriangles,"inputHighOrderTriangles","Triangles where high order points are lying"))
	, inputHighOrderEdgePositions(initData(&inputHighOrderEdgePositions,"inputHighOrderEdgePositions","High order edge points of the mesh "))
    , inputHighOrderTrianglePositions(initData(&inputHighOrderTrianglePositions,"inputHighOrderTrianglePositions","High order triangle points of the mesh "))
    , inputHighOrderTetrahedronPositions(initData(&inputHighOrderTetrahedronPositions,"inputHighOrderTetrahedronPositions","High order tetrahedron points of the mesh "))
{
    addAlias(&d_degree, "order");
}



void HighOrderTetrahedronSetTopologyContainer::init()
{
     d_degree.updateIfDirty(); // make sure m_tetrahedron is up to date
	 d_numberOfTetrahedralPoints.updateIfDirty();
	TetrahedronSetTopologyContainer::init(); // initialize the tetrahedron array
	reinit();
}
size_t HighOrderTetrahedronSetTopologyContainer::getLexicographicIndex(size_t i) const {
	return(lexicographicIndexArray[i]);
}
size_t HighOrderTetrahedronSetTopologyContainer::getHierarchicalIndex(size_t i) const {
	return(hierarchicalIndexArray[i]);
}

size_t HighOrderTetrahedronSetTopologyContainer::getLexicographicIndex(const TetrahedronIndexVector tbi1) const {
	HighOrderDegreeType degree=d_degree.getValue();
	return((size_t)( (degree+1)*(degree+2)*(degree+3)/6-
		 (degree+1-tbi1[0])*(degree+2-tbi1[0])*(degree+3-tbi1[0])/6+
		 (degree+1-tbi1[0])*(degree+2-tbi1[0])/2-
		 (degree+1-tbi1[0]-tbi1[1])*(degree+2-tbi1[0]-tbi1[1])/2+ tbi1[2]));
}
void HighOrderTetrahedronSetTopologyContainer::reinit()
{
	if (d_degree.getValue()>0) {
		// fill the elementMap and the 3 offsetMap in order to get the global index of an element from its Bezier index 
		HighOrderDegreeType degree=d_degree.getValue();
		HighOrderDegreeType i,j,k;
		size_t nbDOFs=(4+6*(degree-1)+2*(degree-1)*(degree-2)+(degree-1)*(degree-2)*(degree-3)/6);
		lexicographicIndexArray.resize(nbDOFs);
		hierarchicalIndexArray.resize(nbDOFs);
		size_t lexicographicalIndex;
		size_t localIndex=0;
		// vertex index
		for (i=0;i<4;++i) {
			TetrahedronIndexVector bti(0,0,0,0);
			bti[i]=degree;
			elementMap.insert(ElementMapType(bti,ElementTetrahedronIndex(i,-1,-1,-1)));
			localIndexMap.insert(OffsetMapType(bti,localIndex));
			lexicographicalIndex=getLexicographicIndex(bti);
			lexicographicIndexArray[tetrahedronIndexArray.size()]=lexicographicalIndex;
			hierarchicalIndexArray[lexicographicalIndex]=tetrahedronIndexArray.size();
			tetrahedronIndexArray.push_back(bti);
			localIndex++;
		}
		// edge index
		if (degree>1) {
			for (i=0;i<6;++i) {

				for (j=1;j<degree;++j) {
					TetrahedronIndexVector bti(0,0,0,0);
					bti[edgesInTetrahedronArray[i][0]]=degree-j;
					bti[edgesInTetrahedronArray[i][1]]=j;
					elementMap.insert(ElementMapType(bti,ElementTetrahedronIndex(-1,i,-1,-1)));
					edgeOffsetMap.insert(OffsetMapType(bti,j-1));
					localIndexMap.insert(OffsetMapType(bti,localIndex));
					lexicographicalIndex=getLexicographicIndex(bti);
					lexicographicIndexArray[tetrahedronIndexArray.size()]=lexicographicalIndex;
					hierarchicalIndexArray[lexicographicalIndex]=tetrahedronIndexArray.size();
					tetrahedronIndexArray.push_back(bti);
					localIndex++;
				}
			}
		}
		// triangle index
		if (degree>2) {
			size_t ind;
			offsetToTriangleIndexArray.clear();
			for (i=0;i<4;++i) {
				for (ind=0,j=1;j<(degree-1);++j) {
					for (k=1;k<(degree-j);++k,++ind) {
						TetrahedronIndexVector bti(0,0,0,0);
						bti[trianglesInTetrahedronArray[i][0]]=j;
						bti[trianglesInTetrahedronArray[i][1]]=k;
						bti[trianglesInTetrahedronArray[i][2]]=degree-j-k;
						if (i==0)
							offsetToTriangleIndexArray.push_back(TriangleIndexVector(j,k,degree-k-j));			
						elementMap.insert(ElementMapType(bti,ElementTetrahedronIndex(-1,-1,i,-1)));
						triangleOffsetMap.insert(OffsetMapType(bti,ind));
						//						std::cerr << "offsetMap["<<(size_t)bti[0]<<' '<<(size_t)bti[1]<<' '<<(size_t)bti[2]<<' '<<(size_t)bti[3]<<" ]= "<<ind<<std::endl;
						localIndexMap.insert(OffsetMapType(bti,localIndex));
						lexicographicalIndex=getLexicographicIndex(bti);
						lexicographicIndexArray[tetrahedronIndexArray.size()]=lexicographicalIndex;
						hierarchicalIndexArray[lexicographicalIndex]=tetrahedronIndexArray.size();
						tetrahedronIndexArray.push_back(bti);
						localIndex++;
					}
				}
			}
		}
		// tetrahedron index
		if (degree>3) {
			size_t ind=0;
			offsetToTetrahedronIndexArray.clear();
			for (i=1;i<(degree-2);++i) {
				for (j=1;j<(degree-i-1);++j) {
					for (k=1;k<(degree-j-i);++k,++ind) {
						TetrahedronIndexVector bti(0,0,0,0);
						bti[0]=i;bti[1]=j;bti[2]=k;
						bti[3]=degree-i-j-k;
						offsetToTetrahedronIndexArray.push_back(bti);
						elementMap.insert(ElementMapType(bti,ElementTetrahedronIndex(-1,-1,-1,0)));
						tetrahedronOffsetMap.insert(OffsetMapType(bti,ind));
						localIndexMap.insert(OffsetMapType(bti,localIndex));
						lexicographicalIndex=getLexicographicIndex(bti);
						lexicographicIndexArray[tetrahedronIndexArray.size()]=lexicographicalIndex;
						hierarchicalIndexArray[lexicographicalIndex]=tetrahedronIndexArray.size();
						tetrahedronIndexArray.push_back(bti);
						localIndex++;
					}
				}
			}
		}
		// initialize the array of weights if necessary

		if (d_weightArray.getValue().empty()){

			SeqWeights &wa=*(d_weightArray.beginEdit());
			wa.resize(this->getNbPoints());
			std::fill(wa.begin(),wa.end(),(SReal)1);
			d_weightArray.endEdit();
		}
		if (d_isRationalSpline.getValue().empty()){
			helper::WriteOnlyAccessor<Data <SeqBools> >  isRationalSpline=d_isRationalSpline;
			isRationalSpline.resize(this->getNbPoints());
			std::fill(isRationalSpline.begin(),isRationalSpline.end(),false);
		}
		// manually creates the edge and triangle structures.
		createTriangleSetArray();
		createEdgeSetArray();
		createEdgesInTetrahedronArray();
		createTrianglesInTetrahedronArray();
	}
	if (inputHighOrderEdgePositions.getValue().size()>0)
		parseInputData();
	if ((d_numberOfTetrahedralPoints.getValue()==0) && (getNumberOfTetrahedra()>0)){
		// compute the number of tetrahedral point if it is not provided
		std::set<size_t> vertexSet;
		size_t i;
		// count the number of vertices involved in the list of tetrahedra
		const sofa::helper::vector<Tetrahedron> &tra=getTetrahedronArray();
		for (i=0;i<tra.size();++i) {
			vertexSet.insert(tra[i][0]);
			vertexSet.insert(tra[i][1]);
			vertexSet.insert(tra[i][2]);
			vertexSet.insert(tra[i][3]);
		}
		d_numberOfTetrahedralPoints.setValue(vertexSet.size());

	}
	// initialize 	weightedDOFArray
	tetrahedronDOFArray.clear(); 
	VecPointID indexArray;
	size_t i;
    for (i=0;i<(size_t)getNbTetrahedra();++i) {
		indexArray.clear();
		getGlobalIndexArrayOfControlPointsInTetrahedron(i,indexArray);
		tetrahedronDOFArray.push_back(indexArray);
	}

//	checkHighOrderTetrahedronTopology();
}
void HighOrderTetrahedronSetTopologyContainer::parseInputData()  {

	HighOrderDegreeType degree=d_degree.getValue();
	size_t i,j,k;
	// count the number of vertices involved in the list of triangles
	const sofa::helper::vector<Tetrahedron> &tta=getTetrahedronArray();
	for (i=0;i<tta.size();++i) {
		for (j=0;j<4;++j) {
			ControlPointLocation cpl(tta[i][j],std::make_pair(POINT,0));
			locationToGlobalIndexMap.insert(std::pair<ControlPointLocation,size_t>(cpl,tta[i][j]));
			globalIndexToLocationMap.insert(std::pair<size_t,ControlPointLocation>(tta[i][j],cpl));
		}

	}


	helper::ReadAccessor<Data<helper::vector< HighOrderEdgePosition > > > hoEdgeArray = this->inputHighOrderEdgePositions;
	helper::ReadAccessor<Data<helper::vector< Edge > > > inputEdgeArray = this->inputEdges;
	for (size_t i = 0; i < hoEdgeArray.size(); ++i)
	{
		size_t pointIndex=hoEdgeArray[i][0];
		size_t edgeIndex=hoEdgeArray[i][1];
		size_t tmp=hoEdgeArray[i][2]+hoEdgeArray[i][3];
		assert((hoEdgeArray[i][2]+hoEdgeArray[i][3])==degree);
		size_t v0=std::min(inputEdgeArray[edgeIndex][0],inputEdgeArray[edgeIndex][1]);
		size_t v1=std::max(inputEdgeArray[edgeIndex][0],inputEdgeArray[edgeIndex][1]);
		int realEdgeIndex=this->getEdgeIndex(v0,v1);
		assert(realEdgeIndex>=0);
		if (realEdgeIndex>=0) {
			size_t offset;
			if (v0==inputEdgeArray[edgeIndex][0]) {
				offset=degree-hoEdgeArray[i][2]-1;
			} else {
				offset=degree-hoEdgeArray[i][3]-1;
			}
			ControlPointLocation cpl((EdgeID)realEdgeIndex,std::make_pair(EDGE,offset));
			locationToGlobalIndexMap.insert(std::pair<ControlPointLocation,size_t>(cpl,pointIndex));
			globalIndexToLocationMap.insert(std::pair<size_t,ControlPointLocation>(pointIndex,cpl));
		}
	}
	
	helper::ReadAccessor<Data<helper::vector< HighOrderTrianglePosition > > > hoTriangleArray = this->inputHighOrderTrianglePositions;
	helper::ReadAccessor<Data<helper::vector< Triangle > > > inputTriangleArray = this->inputTriangles;
	for (size_t i = 0; i < hoTriangleArray.size(); ++i)
	{
		size_t pointIndex= hoTriangleArray[i][0];
		size_t triangleIndex= hoTriangleArray[i][1];
		assert((hoTriangleArray[i][2]+ hoTriangleArray[i][3]+ hoTriangleArray[i][4])==degree);
		Triangle t=inputTriangleArray[triangleIndex];
		// now find a corresponding triangle in the set of triangles
		int realTriangleIndex=this->getTriangleIndex(t[0],t[1],t[2]);
		if (realTriangleIndex>=0) {
			size_t numbering[3];
			Triangle t2=this->getTriangle((TriangleID) realTriangleIndex); 
			TetrahedronIndexVector tvi;
			for (j=0;j<3;++j) {
				for(k=0;t[k]!=t2[j];++k);
				tvi[k]=
				numbering[j]=k;
			}
			size_t offset=(hoTriangleArray[i][numbering[0]+2]-1)*(degree-1)+ hoTriangleArray[i][numbering[1]+2]-1;
			ControlPointLocation cpl((TriangleID) realTriangleIndex,std::make_pair(TRIANGLE,offset));
			locationToGlobalIndexMap.insert(std::pair<ControlPointLocation,size_t>(cpl,pointIndex));
			globalIndexToLocationMap.insert(std::pair<size_t,ControlPointLocation>(pointIndex,cpl));

		}
	}
	
	helper::ReadAccessor<Data<helper::vector< HighOrderTetrahedronPosition > > > hoTetrahedronArray = this->inputHighOrderTetrahedronPositions;
	for (size_t i = 0; i < hoTetrahedronArray.size(); ++i)
	{
		size_t pointIndex= hoTetrahedronArray[i][0];
		size_t tetrahedronIndex= hoTetrahedronArray[i][1];
		assert((hoTetrahedronArray[i][2]+ hoTetrahedronArray[i][3]+ hoTetrahedronArray[i][4]+ hoTetrahedronArray[i][5])==degree);
		Triangle t=this->getTriangle(tetrahedronIndex);
		TetrahedronIndexVector tiv(hoTetrahedronArray[i][2], hoTetrahedronArray[i][3], hoTetrahedronArray[i][4], hoTetrahedronArray[i][5]);
		OffsetMapIterator omi=edgeOffsetMap.find(tiv);
		assert(omi!=edgeOffsetMap.end());
		ControlPointLocation cpl(tetrahedronIndex,std::make_pair(TETRAHEDRON,(*omi).second));
		locationToGlobalIndexMap.insert(std::pair<ControlPointLocation,size_t>(cpl,pointIndex));
		globalIndexToLocationMap.insert(std::pair<size_t,ControlPointLocation>(pointIndex,cpl));	
	}

}
const HighOrderTetrahedronSetTopologyContainer::VecPointID &HighOrderTetrahedronSetTopologyContainer::getGlobalIndexArrayOfControlPoints(const TetraID tetrahedronIndex) const {
	return tetrahedronDOFArray[tetrahedronIndex];
}
HighOrderDegreeType HighOrderTetrahedronSetTopologyContainer::getDegree() const{
	return d_degree.getValue();
}
 size_t HighOrderTetrahedronSetTopologyContainer::getNumberOfTetrahedralPoints() const{
	 return d_numberOfTetrahedralPoints.getValue();
 }
 SReal HighOrderTetrahedronSetTopologyContainer::getWeight(int i) const {
		 return(d_weightArray.getValue()[i]);
 }
 bool HighOrderTetrahedronSetTopologyContainer::isRationalSpline(int i) const {
	 return(d_isRationalSpline.getValue()[i]);
 }
 const HighOrderTetrahedronSetTopologyContainer::SeqWeights & HighOrderTetrahedronSetTopologyContainer::getWeightArray() const {
	 return(d_weightArray.getValue());
 }
size_t HighOrderTetrahedronSetTopologyContainer::getGlobalIndexOfControlPoint(const TetraID tetrahedronIndex,
	const TetrahedronIndexVector id) 
{		
	Tetrahedron tet=getTetrahedron(tetrahedronIndex);
	HighOrderDegreeType degree=d_degree.getValue();
	ElementMapIterator emi=elementMap.find(id);
	if (emi!=elementMap.end()) {
		ElementTetrahedronIndex ei=(*emi).second;
		if (ei[0]!= -1) {
			// point is a vertex of the tetrahedral mesh
			return (tet[ei[0]]);
		} else if (ei[1]!= -1) {
			// point is on an edge of the tetrahedral mesh

			// the points on edges are stored after the tetrahedron vertices
			// there are (degree-1) points store on each edge
			// eit[ei[1]] = id of edge where the point is located
			// ei[1] = the local index (<6) of the edge where the point is located 
			EdgesInTetrahedron eit=getEdgesInTetrahedron(tetrahedronIndex);
			Edge e=getEdge(eit[ei[1]]);
			// test if the edge is along the right direction
			OffsetMapIterator omi=edgeOffsetMap.find(id);
			if (locationToGlobalIndexMap.empty()) {
				if (e[0]==tet[edgesInTetrahedronArray[ei[1]][0]]) {
					return(getNumberOfTetrahedralPoints()+eit[ei[1]]*(degree-1)+(*omi).second);
				} else {
					// use the other direction
					return(getNumberOfTetrahedralPoints()+eit[ei[1]]*(degree-1)+degree-2-(*omi).second);
				}
			} else {
				ControlPointLocation cpl(eit[ei[1]],std::make_pair(EDGE,(*omi).second));
				if (e[0]==tet[edgesInTetrahedronArray[ei[1]][0]]) {
					cpl=ControlPointLocation(eit[ei[1]],std::make_pair(EDGE,degree-2-(*omi).second));
				}
				assert(locationToGlobalIndexMap.find(cpl)!=locationToGlobalIndexMap.end());
				return(locationToGlobalIndexMap[cpl]);
			}
		} else if (ei[2]!= -1) {
			// point is on an edge of the tetrahedral mesh
			TrianglesInTetrahedron tit=getTrianglesInTetrahedron(tetrahedronIndex);
			// the points on triangles are stored after the Bezier points on edges
			// there are (degree-1)*(degree-2)/2 points store on each triangle
			// eit[ei[2]] = id of triangle where the point is located
			// ei[2] = the local index (<4) of the triangle where the point is located 
			Triangle tr=getTriangle(tit[ei[2]]);
			Triangle indexTriangle;
			size_t k,i;
			i=ei[2];
			for (k=0;(tr[0]!=tet[trianglesInTetrahedronArray[i][k]]);++k);
			//				indexTriangle[0]=k;
			indexTriangle[k]=0;
			if (tr[1]==tet[trianglesInTetrahedronArray[i][(k+1)%3]]) {
				indexTriangle[(k+1)%3]=1;
				indexTriangle[(k+2)%3]=2;
			} else {
				indexTriangle[(k+2)%3]=1;
				indexTriangle[(k+1)%3]=2;
			}
			TetrahedronIndexVector bti(0,0,0,0);
			bti[trianglesInTetrahedronArray[i][indexTriangle[0]]]=id[trianglesInTetrahedronArray[i][0]];
			bti[trianglesInTetrahedronArray[i][indexTriangle[1]]]=id[trianglesInTetrahedronArray[i][1]];
			bti[trianglesInTetrahedronArray[i][indexTriangle[2]]]=id[trianglesInTetrahedronArray[i][2]];
			OffsetMapIterator omi=triangleOffsetMap.find(bti);
			if (locationToGlobalIndexMap.empty()) {
				return(getNumberOfTetrahedralPoints()+getNumberOfEdges()*(degree-1)+tit[i]*(degree-1)*(degree-2)/2+(*omi).second);
			} else {
				ControlPointLocation cpl(tit[i],std::make_pair(TRIANGLE,(*omi).second));
				assert(locationToGlobalIndexMap.find(cpl)!=locationToGlobalIndexMap.end());
				return(locationToGlobalIndexMap[cpl]);
			}
		} else if (ei[3]!= -1) {
			// the points on edges are stored after the tetrahedron vertices
			// there are (degree-1)*(degree-2)/2 points store on each edge
			// eit[ei[1]] = id of edge where the point is located
			// ei[1] = the local index (<6) of the edge where the point is located 
			OffsetMapIterator omi=tetrahedronOffsetMap.find(id);
			if (locationToGlobalIndexMap.empty()) {
				return(getNumberOfTetrahedralPoints()+getNumberOfEdges()*(degree-1)+getNumberOfTriangles()*(degree-1)*(degree-2)/2+tetrahedronIndex*(degree-1)*(degree-2)*(degree-3)/6+(*omi).second);
			} else {
				ControlPointLocation cpl(ei[3],std::make_pair(TETRAHEDRON,(*omi).second));
				assert(locationToGlobalIndexMap.find(cpl)!=locationToGlobalIndexMap.end());
				return(locationToGlobalIndexMap[cpl]);
			}
		}
	} else {
#ifndef NDEBUG
		sout << "Error. [HighOrderTetrahedronSetTopologyContainer::getGlobalIndexOfControlPoint] Bezier Index "<< id <<" has not been recognized to be valid" << sendl;
#endif
        return (0);
    }
    return (0);
}

TetrahedronIndexVector HighOrderTetrahedronSetTopologyContainer::getTetrahedronIndex(const size_t localIndex) const
{
	
	if (localIndex<tetrahedronIndexArray.size()) {
		return tetrahedronIndexArray[localIndex];
	} else {
#ifndef NDEBUG
		sout << "Error. [HighOrderTetrahedronSetTopologyContainer::getTetrahedronIndexVector] Index "<< localIndex <<" is greater than the number "<< tetrahedronIndexArray.size() <<" of control points." << sendl;
#endif
		TetrahedronIndexVector id;
		return (id);
	}
}
size_t HighOrderTetrahedronSetTopologyContainer::getLocalIndexFromTetrahedronIndex(const TetrahedronIndexVector id) const {
	OffsetMapConstIterator omi=localIndexMap.find(id);
	if (omi==localIndexMap.end())
	{
#ifndef NDEBUG
		sout << "Error. [HighOrderTetrahedronSetTopologyContainer::getLocalIndexFromTetrahedronIndexVector] Tetrahedron Bezier Index "<< id  <<" is out of range." << sendl;
#endif
		return(0);
	} else {
		return ((*omi).second);
	}
}
sofa::helper::vector<TetrahedronIndexVector> HighOrderTetrahedronSetTopologyContainer::getTetrahedronIndexArray() const
{
	return (tetrahedronIndexArray);
}
sofa::helper::vector<TetrahedronIndexVector> HighOrderTetrahedronSetTopologyContainer::getTetrahedronIndexArrayOfGivenDegree(const HighOrderDegreeType deg) const
{
	// vertex index
	size_t i,j,k;
	sofa::helper::vector<TetrahedronIndexVector> tbiArray;
	for (i=0;i<4;++i) {
		TetrahedronIndexVector bti(0,0,0,0);
		bti[i]=deg;
		tbiArray.push_back(bti);
	}
	// edge index
	if (deg>1) {
		for (i=0;i<6;++i) {
			for (j=1;j<deg;++j) {
				TetrahedronIndexVector bti(0,0,0,0);
				bti[edgesInTetrahedronArray[i][0]]=deg-j;
				bti[edgesInTetrahedronArray[i][1]]=j;
				tbiArray.push_back(bti);
			}
		}
	}
	// triangle index
	if (deg>2) {;
		for (i=0;i<4;++i) {
			for (j=1;j<(size_t)(deg-1);++j) {
				for (k=1;k<(deg-j);++k) {
					TetrahedronIndexVector bti(0,0,0,0);
					bti[trianglesInTetrahedronArray[i][0]]=j;
					bti[trianglesInTetrahedronArray[i][1]]=k;
					bti[trianglesInTetrahedronArray[i][2]]=deg-j-k;
					tbiArray.push_back(bti);
				}
			}
		}
	}
	// tetrahedron index
	if (deg>3) {
		for (i=1;i<(size_t)(deg-2);++i) {
			for (j=1;j<(size_t)(deg-i-1);++j) {
				for (k=1;k<(deg-j-i);++k) {
					TetrahedronIndexVector bti(0,0,0,0);
					bti[0]=i;bti[1]=j;bti[2]=k;
					bti[3]=deg-i-j-k;
					tbiArray.push_back(bti);
				}
			}
		}
	}
	return(tbiArray);
}
sofa::helper::vector<HighOrderTetrahedronSetTopologyContainer::LocalTetrahedronIndex> HighOrderTetrahedronSetTopologyContainer::getMapOfTetrahedronIndexArrayFromInferiorDegree() const
{
	HighOrderDegreeType degree=d_degree.getValue();
	sofa::helper::vector<TetrahedronIndexVector> tbiDerivArray=getTetrahedronIndexArrayOfGivenDegree(degree-1);
	sofa::helper::vector<TetrahedronIndexVector> tbiLinearArray=getTetrahedronIndexArrayOfGivenDegree(1);
	TetrahedronIndexVector tbi;
	sofa::helper::vector<LocalTetrahedronIndex> correspondanceArray;
//	correspondanceArray.resize(tbiDerivArray.size());
	size_t i,j;
	for (i=0;i<tbiDerivArray.size();++i) {
		LocalTetrahedronIndex correspondance;
		for (j=0;j<4;++j) {
			tbi=tbiDerivArray[i]+tbiLinearArray[j];
			correspondance[j]=getLocalIndexFromTetrahedronIndex(tbi);
		}
		correspondanceArray.push_back(correspondance);
	}
	return(correspondanceArray);
}
void HighOrderTetrahedronSetTopologyContainer::getGlobalIndexArrayOfControlPointsInTetrahedron(const TetraID tetrahedronIndex, VecPointID & indexArray)
{
	Tetrahedron tet=getTetrahedron(tetrahedronIndex);
	indexArray.clear();
	// vertex index
	size_t i,j,k;
	for (i=0;i<4;++i)
		indexArray.push_back(tet[i]);

	size_t offset;
	// edge index
	HighOrderDegreeType degree=d_degree.getValue();
	if (degree>1) {
		EdgesInTetrahedron eit=getEdgesInTetrahedron(tetrahedronIndex);
		for (i=0;i<6;++i) {
			Edge e=getEdge(eit[i]);
			if (locationToGlobalIndexMap.empty()) {
				offset=getNumberOfTetrahedralPoints()+eit[i]*(degree-1);
				// check the order of the edge to be consistent with the tetrahedron
				if (e[0]==tet[edgesInTetrahedronArray[i][0]]) {
					for (j=0;j<(size_t)(degree-1);++j) {
						indexArray.push_back(offset+j);
					}
				} else {
					int jj;
					for (jj=degree-2;jj>=0;--jj) {
						indexArray.push_back(offset+jj);
					}
				}
			} else {
				ControlPointLocation cpl;
				
				if (e[0]!=tet[edgesInTetrahedronArray[i][0]]) {
					for (j=0;j<(size_t)(degree-1);++j) {
						cpl=ControlPointLocation(eit[i],std::make_pair(EDGE,j));
						assert(locationToGlobalIndexMap.find(cpl)!=locationToGlobalIndexMap.end());
						indexArray.push_back(locationToGlobalIndexMap[cpl]);
					}
				}
				else {
					int jj;
					for (jj=degree-2;jj>=0;--jj) {
						cpl=ControlPointLocation(eit[i],std::make_pair(EDGE,(size_t)jj));
						assert(locationToGlobalIndexMap.find(cpl)!=locationToGlobalIndexMap.end());
						indexArray.push_back(locationToGlobalIndexMap[cpl]);
					}
				}

			}
		}
	}
	// triangle index
	if (degree>2) {
		TrianglesInTetrahedron tit=getTrianglesInTetrahedron(tetrahedronIndex);
		size_t pointsPerTriangle=(degree-1)*(degree-2)/2;
		for (i=0;i<4;++i) {
			offset=getNumberOfTetrahedralPoints()+getNumberOfEdges()*(degree-1)+tit[i]*pointsPerTriangle;
			Triangle tr=getTriangle(tit[i]);
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
			for (j=1;j<(size_t)(degree-1);++j) {
				for (k=1;k<(degree-j);++k) {
					TetrahedronIndexVector bti(0,0,0,0);
					bti[trianglesInTetrahedronArray[i][indexTriangle[0]]]=j;
					bti[trianglesInTetrahedronArray[i][indexTriangle[1]]]=k;
					bti[trianglesInTetrahedronArray[i][indexTriangle[2]]]=degree-j-k;
					OffsetMapIterator omi=triangleOffsetMap.find(bti);
					if (locationToGlobalIndexMap.empty()) {
						indexArray.push_back(offset+(*omi).second);
					} else {
						ControlPointLocation cpl(tit[i],std::make_pair(TRIANGLE,(*omi).second));
						assert(locationToGlobalIndexMap.find(cpl)!=locationToGlobalIndexMap.end());
						indexArray.push_back(locationToGlobalIndexMap[cpl]);
					}
				}
			}
		}
	}
	

	// tetrahedron index
	if (degree>3) {
		size_t pointsPerTetrahedron=(degree-1)*(degree-2)*(degree-3)/6;
		offset=getNumberOfTetrahedralPoints()+getNumberOfEdges()*(degree-1)+getNumberOfTriangles()*(degree-1)*(degree-2)/2+tetrahedronIndex*pointsPerTetrahedron;
		size_t rank=0;
		for (i=0;i<(size_t)(degree-3);++i) {
			for (j=0;j<(degree-i-3);++j) {
				for (k=0;k<(degree-j-i-3);++k) {
					if (locationToGlobalIndexMap.empty()) {
						indexArray.push_back(offset+rank);
					} else {
						ControlPointLocation cpl(tetrahedronIndex,std::make_pair(TETRAHEDRON,rank));
						assert(locationToGlobalIndexMap.find(cpl)!=locationToGlobalIndexMap.end());
						indexArray.push_back(locationToGlobalIndexMap[cpl]);
					}
					rank++;
				}
			}
		}
	}

	}
void HighOrderTetrahedronSetTopologyContainer::getLocationFromGlobalIndex(const size_t globalIndex, HighOrderTetrahedronPointLocation &location, 
	size_t &elementIndex, size_t &elementOffset)
{
	size_t gi=globalIndex;
	if (globalIndexToLocationMap.empty()) {
		if (gi<getNumberOfTetrahedralPoints()) {
			location=POINT;
			elementIndex=gi;
			elementOffset=0;
		} else {
			gi-=getNumberOfTetrahedralPoints();
			HighOrderDegreeType degree=d_degree.getValue();
			if (gi<(getNumberOfEdges()*(degree-1))) {
				location=EDGE;
				elementIndex=gi/(degree-1);
				elementOffset=gi%(degree-1);
			} else {
				gi-=getNumberOfEdges()*(degree-1);
				size_t pointsPerTriangle=(degree-1)*(degree-2)/2;
				if (gi<(getNumberOfTriangles()*pointsPerTriangle)) {
					location=TRIANGLE;
					elementIndex=gi/pointsPerTriangle;
					elementOffset=gi%pointsPerTriangle;
				} else {
					gi-=getNumberOfTriangles()*pointsPerTriangle;
					size_t pointsPerTetrahedron=(degree-1)*(degree-2)*(degree-3)/6;
					if (gi<(getNumberOfTetrahedra()*pointsPerTetrahedron)) {
						location=TETRAHEDRON;
						elementIndex=gi/pointsPerTetrahedron;
						elementOffset=gi%pointsPerTetrahedron;
					}  else {
#ifndef NDEBUG
						sout << "Error. [HighOrderTetrahedronSetTopologyContainer::getLocationFromGlobalIndex] Global Index "<< globalIndex <<" exceed the number of Bezier Points" << sendl;
#endif
					}
				}
			}
		}
	} else {
		std::map<size_t,ControlPointLocation>::iterator itcpl=globalIndexToLocationMap.find(globalIndex);
#ifndef NDEBUG
		if (itcpl==globalIndexToLocationMap.end())
			sout << "Error. [HighOrderTetrahedronSetTopologyContainer::getLocationFromGlobalIndex] Global Index "<< globalIndex <<" is not in the map globalIndexToLocationMap " << sendl;
		assert(itcpl!=globalIndexToLocationMap.end());
#endif
		location=((*itcpl).second).second.first;
		elementIndex=((*itcpl).second).first;
		elementOffset=((*itcpl).second).second.second;
	}
}
void HighOrderTetrahedronSetTopologyContainer::getEdgeIndexFromEdgeOffset(size_t offset, EdgeIndexVector &ebi){
	assert(offset<d_degree.getValue());
	ebi[0]=offset+1;
	ebi[1]=d_degree.getValue()-offset-1;
}
void HighOrderTetrahedronSetTopologyContainer::getTriangleIndexFromTriangleOffset(size_t offset, TriangleIndexVector &tbi){
    assert(offset<(size_t)((d_degree.getValue()-1)*(d_degree.getValue()-2)/2));
	tbi=offsetToTriangleIndexArray[offset];
}
void HighOrderTetrahedronSetTopologyContainer::getTetrahedronIndexFromTetrahedronOffset(size_t offset, TetrahedronIndexVector &tbi){
    assert(offset<(size_t)((d_degree.getValue()-1)*(d_degree.getValue()-2)*(d_degree.getValue()-3)/6));
	tbi=offsetToTetrahedronIndexArray[offset];
}
bool HighOrderTetrahedronSetTopologyContainer::checkHighOrderTetrahedronTopology()
{
	size_t nTetras,elem;
	HighOrderDegreeType degree=d_degree.getValue();
	// check the total number of vertices.
    assert((size_t)getNbPoints()==(getNumberOfTetrahedralPoints()+getNumberOfEdges()*(degree-1)+getNumberOfTriangles()*(degree-1)*(degree-2)/2+getNumberOfTetrahedra()*(degree-1)*(degree-2)*(degree-3)/6));
	sofa::helper::vector<TetrahedronIndexVector> tbiArray=getTetrahedronIndexArray();
	VecPointID indexArray;
	HighOrderTetrahedronPointLocation location; 
    size_t elementIndex, elementOffset/*,localIndex*/;
	for (nTetras=0;nTetras<getNumberOfTetrahedra();++nTetras) {
		indexArray.clear();
		getGlobalIndexArrayOfControlPointsInTetrahedron(nTetras,indexArray);
		// check the number of control points per tetrahedron is correct
        assert(indexArray.size()==(size_t)((4+6*(degree-1)+2*(degree-1)*(degree-2)+(degree-1)*(degree-2)*(degree-3)/6)));
		for(elem=0;elem<indexArray.size();++elem) {
			size_t globalIndex=getGlobalIndexOfControlPoint(nTetras,tbiArray[elem]);
			// check that getGlobalIndexOfControlPoint and getGlobalIndexArrayOfBezierPointsInTetrahedron give the same answer
			assert(globalIndex==indexArray[elem]);
#ifndef NDEBUG
            TetrahedronIndexVector tbi=getTetrahedronIndex(elem);
#endif
			assert(elem==getLocalIndexFromTetrahedronIndex(tbi));
			// check that getTetrahedronIndexVector is consistant with getTetrahedronIndexVectorArray
			assert(tbiArray[elem][0]==tbi[0]);
			assert(tbiArray[elem][1]==tbi[1]);
			assert(tbiArray[elem][2]==tbi[2]);
			assert(tbiArray[elem][3]==tbi[3]);
			// check that getLocationFromGlobalIndex is consistent with 
			getLocationFromGlobalIndex(globalIndex,location,elementIndex,elementOffset);
			if (elem<4) {
				assert(location==POINT);
				assert(elementIndex==getTetrahedron(nTetras)[elem]);
				assert(elementOffset==0);
			}
			else if (elem<(size_t)(4+6*(degree-1))){
				assert(location==EDGE);
				assert(elementIndex==getEdgesInTetrahedron(nTetras)[(elem-4)/(degree-1)]);
			}
			else if (elem<(size_t)(4+6*(degree-1)+2*(degree-1)*(degree-2))){
                assert(location==TRIANGLE);
#ifndef NDEBUG
                size_t nbPointPerEdge=(degree-1)*(degree-2)/2;
                size_t val=(elem-4-6*(degree-1))/(nbPointPerEdge);
#endif
				assert(elementIndex==getTrianglesInTetrahedron(nTetras)[val]);
			}
		}

	}
	if (locationToGlobalIndexMap.size()>0) {
		// check consistency between both maps
		assert(locationToGlobalIndexMap.size()==globalIndexToLocationMap.size());
		std::map<ControlPointLocation,size_t>::iterator itcpl;
		std::map<size_t,ControlPointLocation>::iterator itgi;
		for (itcpl=locationToGlobalIndexMap.begin();itcpl!=locationToGlobalIndexMap.end();++itcpl) {
			itgi=globalIndexToLocationMap.find(itcpl->second);
			assert(itgi!=globalIndexToLocationMap.end());
			assert((*itgi).second==(*itcpl).first);
		}
	}
	return( true);
}

} // namespace topology

} // namespace component

} // namespace sofa
