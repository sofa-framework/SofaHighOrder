#include "HighOrderTriangleSetTopologyContainer.h"
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

SOFA_DECL_CLASS(HighOrderTriangleSetTopologyContainer)
int HighOrderTriangleSetTopologyContainerClass = core::RegisterObject("Bezier Triangle set topology container")
        .add< HighOrderTriangleSetTopologyContainer >()
        ;


typedef std::pair<TriangleIndexVector,HighOrderTriangleSetTopologyContainer::ElementTriangleIndex> ElementMapType;
typedef std::map<TriangleIndexVector,HighOrderTriangleSetTopologyContainer::ElementTriangleIndex>::iterator ElementMapIterator;

typedef std::pair<TriangleIndexVector,size_t> OffsetMapType;
typedef std::map<TriangleIndexVector,size_t>::iterator OffsetMapIterator;
typedef std::map<TriangleIndexVector,size_t>::const_iterator OffsetMapConstIterator;

HighOrderTriangleSetTopologyContainer::HighOrderTriangleSetTopologyContainer()
    : TriangleSetTopologyContainer()
    , d_degree(initData(&d_degree, (size_t)0,"degree", "Degree of Bezier Tetrahedra"))
    , d_numberOfTriangularPoints(initData(&d_numberOfTriangularPoints, (size_t) 0,"NbTriangularVertices", "Number of Triangular Vertices"))
    , d_isRationalSpline(initData(&d_isRationalSpline, SeqBools(),"isRational", "If a bezier triangle  is rational or integral"))
	, d_weightArray(initData(&d_weightArray, SeqWeights(),"weights", "Array of weights for rational bezier triangles"))
    , inputEdges(initData(&inputEdges,"inputHighOrderEdges","Edges where high order points are lying"))
	, inputHighOrderEdgePositions(initData(&inputHighOrderEdgePositions,"inputHighOrderEdgePositions","High order edge points of the mesh loaded"))
    , inputHighOrderTrianglePositions(initData(&inputHighOrderTrianglePositions,"inputHighOrderTrianglePositions","High order triangle points of the mesh loaded"))
{
    addAlias(&d_degree, "order");
}



void HighOrderTriangleSetTopologyContainer::init()
{
     d_degree.updateIfDirty(); // make sure m_Triangle is up to date
	 d_numberOfTriangularPoints.updateIfDirty();
	TriangleSetTopologyContainer::init(); // initialize the Triangle array
	reinit();
}
void HighOrderTriangleSetTopologyContainer::reinit()
{
 	if (d_degree.getValue()>0) {
		// clear previous entries if it exists
		elementMap.clear();
		localIndexMap.clear();
		bezierIndexArray.clear();
		edgeOffsetMap.clear();
		triangleOffsetMap.clear();

		// fill the elementMap and the 3 offsetMap in order to get the global index of an element from its Bezier index 
		HighOrderDegreeType degree=d_degree.getValue();
		HighOrderDegreeType i,j;
		size_t nbDOFs=3+3*(degree-1)+(degree-1)*(degree-2)/2;
		lexicographicIndexArray.resize(nbDOFs);
		hierarchicalIndexArray.resize(nbDOFs);
		size_t lexicographicalIndex;

		size_t localIndex=0;
		// vertex index
		for (i=0;i<3;++i) {
			TriangleIndexVector bti(0,0,0);
			bti[i]=degree;
			elementMap.insert(ElementMapType(bti,ElementTriangleIndex(i,-1,-1)));
			localIndexMap.insert(OffsetMapType(bti,localIndex));
			lexicographicalIndex=getLexicographicIndex(bti);
			lexicographicIndexArray[bezierIndexArray.size()]=lexicographicalIndex;
			hierarchicalIndexArray[lexicographicalIndex]=bezierIndexArray.size();
			bezierIndexArray.push_back(bti);
			localIndex++;
		}
		// edge index
		if (degree>1) {
			for (i=0;i<3;++i) {

				for (j=1;j<degree;++j) {
					TriangleIndexVector bti(0,0,0);
					bti[(i+1)%3]=degree-j;
					bti[(i+2)%3]=j;
					elementMap.insert(ElementMapType(bti,ElementTriangleIndex(-1,i,-1)));
					edgeOffsetMap.insert(OffsetMapType(bti,j-1));
					localIndexMap.insert(OffsetMapType(bti,localIndex));
					lexicographicalIndex=getLexicographicIndex(bti);
					lexicographicIndexArray[bezierIndexArray.size()]=lexicographicalIndex;
					hierarchicalIndexArray[lexicographicalIndex]=bezierIndexArray.size();
					bezierIndexArray.push_back(bti);
					localIndex++;
				}
			}
		}
		// triangle index
		if (degree>2) {
			offsetToTriangleIndexVectorArray.clear();
			size_t ind=0;
			for (i=1;i<(degree-1);++i) {
				for (j=1;j<(degree-i);++j,++ind) {
					TriangleIndexVector bti(0,0,0);
					bti[0]=i;bti[1]=j;
					bti[2]=degree-i-j;
					offsetToTriangleIndexVectorArray.push_back(bti);
					elementMap.insert(ElementMapType(bti,ElementTriangleIndex(-1,-1,0)));
					triangleOffsetMap.insert(OffsetMapType(bti,ind));
					//						std::cerr << "offsetMap["<<(size_t)bti[0]<<' '<<(size_t)bti[1]<<' '<<(size_t)bti[2]<<' '<<(size_t)bti[3]<<" ]= "<<ind<<std::endl;
					localIndexMap.insert(OffsetMapType(bti,localIndex));
					lexicographicalIndex=getLexicographicIndex(bti);
					lexicographicIndexArray[bezierIndexArray.size()]=lexicographicalIndex;
					hierarchicalIndexArray[lexicographicalIndex]=bezierIndexArray.size();
					bezierIndexArray.push_back(bti);
					localIndex++;
				}

			}

		}
		
		// manually creates the edge and triangle structures.
        createEdgeSetArray();
        createEdgesInTriangleArray();
	}
	if (inputHighOrderEdgePositions.getValue().size()>0)
		parseInputData();

	if ((d_numberOfTriangularPoints.getValue()==0) && (getNumberOfTriangles()>0)){
		// compute the number of triangular point if it is not provided
		std::set<size_t> vertexSet;
		size_t i;
		// count the number of vertices involved in the list of triangles
		const sofa::helper::vector<Triangle> &tra=getTriangleArray();
		for (i=0;i<tra.size();++i) {
			vertexSet.insert(tra[i][0]);
			vertexSet.insert(tra[i][1]);
			vertexSet.insert(tra[i][2]);
		}
		d_numberOfTriangularPoints.setValue(vertexSet.size());

	}
	// initialize 	triangleDOFArray
	triangleDOFArray.clear(); 
	VecPointID indexArray;
	size_t i;
	for (i=0;i<(size_t)getNbTriangles();++i) {
		indexArray.clear();
		getGlobalIndexArrayOfControlPointsInTriangle(i,indexArray);
		triangleDOFArray.push_back(indexArray);
	}
    // initialize the array of weights if necessary
    if ((d_isRationalSpline.getValue().empty()) && (getNumberOfTriangles()>0)) {
        helper::WriteOnlyAccessor<Data <SeqBools> >  isRationalSpline = d_isRationalSpline;
        isRationalSpline.resize(this->getNumberOfTriangles());
        if (d_weightArray.getValue().empty()) {
            SeqWeights &wa = *(d_weightArray.beginEdit());
            wa.resize(this->getNbPoints());
            std::fill(wa.begin(), wa.end(), (SReal)1);
            d_weightArray.endEdit();
            // if no weights are provided then all tetrahedra are integral
            std::fill(isRationalSpline.begin(), isRationalSpline.end(), false);
        }
        else {
            size_t j;
            const SeqWeights &wa = d_weightArray.getValue();
            // if weights are provided triangles having one weight different from 1 are label as rational.
            const sofa::helper::vector<Triangle> &tta = getTriangleArray();
            for (i = 0; i < tta.size(); ++i) {
                isRationalSpline[i] = false;
                getGlobalIndexArrayOfControlPointsInTriangle(i, indexArray);
                for (j = 0; j < indexArray.size(); ++j) {
                    if (wa[indexArray[j]] != (SReal)1.0)
                        isRationalSpline[i] = true;
                }
            }
        }
    }

//	checkHighOrderTriangleTopology();
}
void HighOrderTriangleSetTopologyContainer::parseInputData()  {
	std::set<size_t> vertexSet;
	HighOrderDegreeType degree=d_degree.getValue();
	size_t i,j;
	// count the number of vertices involved in the list of triangles
	const sofa::helper::vector<Triangle> &tra=getTriangleArray();
	for (i=0;i<tra.size();++i) {
		for (j=0;j<3;++j) {

			TriangleIndexVector tiv(0,0,0);
			tiv[j]=degree;
			locationToGlobalIndexMap.insert(std::pair<HighOrderTriangleSetTopologyContainer::ControlPointLocation,size_t>(HighOrderTriangleSetTopologyContainer::ControlPointLocation(i,tiv),tra[i][j]));
			globalIndexToLocationMap.insert(std::pair<size_t,HighOrderTriangleSetTopologyContainer::ControlPointLocation>(tra[i][j],HighOrderTriangleSetTopologyContainer::ControlPointLocation(i,tiv)));

		}

	}


	helper::ReadAccessor<Data<helper::vector< HighOrderEdgePosition > > > hoEdgeArray = this->inputHighOrderEdgePositions;
	helper::ReadAccessor<Data<helper::vector< Edge > > > inputEdgeArray = this->inputEdges;
	for (size_t i = 0; i < hoEdgeArray.size(); ++i)
	{
		size_t pointIndex=hoEdgeArray[i][0];
		size_t edgeIndex=hoEdgeArray[i][1];
		assert((hoEdgeArray[i][2]+hoEdgeArray[i][3])==degree);
		size_t v0=std::min(inputEdgeArray[edgeIndex][0],inputEdgeArray[edgeIndex][1]);
		size_t v1=std::max(inputEdgeArray[edgeIndex][0],inputEdgeArray[edgeIndex][1]);
		int realEdgeIndex=this->getEdgeIndex(v0,v1);
		assert(realEdgeIndex>=0);
		if (realEdgeIndex>=0) {
			const TrianglesAroundEdge &tae=this->getTrianglesAroundEdge((EdgeID)realEdgeIndex);
			assert(tae.size()>0);
			for (j=0;j<tae.size();++j) {
				size_t triangleIndex=tae[j];
				int edgeTriangleIndex=this->getEdgeIndexInTriangle(this->getEdgesInTriangle(triangleIndex),(EdgeID)realEdgeIndex);
				assert(edgeTriangleIndex>=0);
				Triangle t=this->getTriangle(triangleIndex);
				TriangleIndexVector tiv(0,0,0);
				if (t[(edgeTriangleIndex+1)%3]==inputEdgeArray[edgeIndex][0]) {
					tiv[(edgeTriangleIndex+1)%3]=hoEdgeArray[i][2];
					tiv[(edgeTriangleIndex+2)%3]=hoEdgeArray[i][3];
					assert(t[(edgeTriangleIndex+2)%3]==inputEdgeArray[edgeIndex][1]);
				} else {
					tiv[(edgeTriangleIndex+2)%3]=hoEdgeArray[i][2];
					tiv[(edgeTriangleIndex+1)%3]=hoEdgeArray[i][3];
					assert(t[(edgeTriangleIndex+1)%3]==inputEdgeArray[edgeIndex][1]);
					assert(t[(edgeTriangleIndex+2)%3]==inputEdgeArray[edgeIndex][0]);
				}
				locationToGlobalIndexMap.insert(std::pair<HighOrderTriangleSetTopologyContainer::ControlPointLocation,size_t>(HighOrderTriangleSetTopologyContainer::ControlPointLocation(triangleIndex,tiv),pointIndex));
				globalIndexToLocationMap.insert(std::pair<size_t,HighOrderTriangleSetTopologyContainer::ControlPointLocation>(pointIndex,HighOrderTriangleSetTopologyContainer::ControlPointLocation(triangleIndex,tiv)));
			}
		}
	}
	helper::ReadAccessor<Data<helper::vector< HighOrderTrianglePosition > > > hoTriangleArray = this->inputHighOrderTrianglePositions;
	for (size_t i = 0; i < hoTriangleArray.size(); ++i)
	{
		size_t pointIndex= hoTriangleArray[i][0];
		size_t triangleIndex= hoTriangleArray[i][1];
		assert((hoTriangleArray[i][2]+ hoTriangleArray[i][3]+ hoTriangleArray[i][4])==degree);
		Triangle t=this->getTriangle(triangleIndex);
		TriangleIndexVector tiv(hoTriangleArray[i][2], hoTriangleArray[i][3], hoTriangleArray[i][4]);
		locationToGlobalIndexMap.insert(std::pair<HighOrderTriangleSetTopologyContainer::ControlPointLocation,size_t>(HighOrderTriangleSetTopologyContainer::ControlPointLocation(triangleIndex,tiv),pointIndex));
		globalIndexToLocationMap.insert(std::pair<size_t,HighOrderTriangleSetTopologyContainer::ControlPointLocation>(pointIndex,HighOrderTriangleSetTopologyContainer::ControlPointLocation(triangleIndex,tiv)));
	}

}
SReal HighOrderTriangleSetTopologyContainer::getWeight(int i) const {
	return(d_weightArray.getValue()[i]);
}
bool HighOrderTriangleSetTopologyContainer::isRationalSpline(int i) const {
	 return(d_isRationalSpline.getValue()[i]);
}
const HighOrderTriangleSetTopologyContainer::SeqWeights & HighOrderTriangleSetTopologyContainer::getWeightArray() const {
	return(d_weightArray.getValue());
}

HighOrderDegreeType HighOrderTriangleSetTopologyContainer::getDegree() const{
	return d_degree.getValue();
}

size_t HighOrderTriangleSetTopologyContainer::getNumberOfTriangularPoints() const{
     return d_numberOfTriangularPoints.getValue();
}

size_t HighOrderTriangleSetTopologyContainer::getGlobalIndexOfControlPoint(const TetraID triangleIndex,
     const TriangleIndexVector id) {

         if (locationToGlobalIndexMap.empty()) {
             Triangle tr=getTriangle(triangleIndex);
             HighOrderDegreeType degree=d_degree.getValue();
             ElementMapIterator emi=elementMap.find(id);
             if (emi!=elementMap.end()) {
                 ElementTriangleIndex ei=(*emi).second;
                 if (ei[0]!= -1) {
                     // point is a vertex of the triangular mesh
                     return (tr[ei[0]]);
                 } else if (ei[1]!= -1) {
                     // point is on an edge of the triangular mesh
                     // the points on edges are stored after the Triangle vertices
                     // there are (degree-1) points store on each edge
                     // eit[ei[1]] = id of edge where the point is located
                     // ei[1] = the local index (<3) of the edge where the point is located
                     EdgesInTriangle eit=getEdgesInTriangle(triangleIndex);
                     Edge e=getEdge(eit[ei[1]]);
                     // test if the edge is along the right direction
                     OffsetMapIterator omi=edgeOffsetMap.find(id);
                     if (e[0]==tr[(ei[1]+1)%3]) {
                         return(getNumberOfTriangularPoints()+eit[ei[1]]*(degree-1)+(*omi).second);
                     } else {
                         // use the other direction
                         return(getNumberOfTriangularPoints()+eit[ei[1]]*(degree-1)+degree-2-(*omi).second);
                     }

                 } else if (ei[2]!= -1) {
                     // the points on edges are stored after the tetrahedron vertices
                     // there are (degree-1)*(degree-2)/2 points store on each edge
                     // eit[ei[1]] = id of edge where the point is located
                     // ei[1] = the local index (<6) of the edge where the point is located
                     OffsetMapIterator omi=triangleOffsetMap.find(id);
                     return(getNumberOfTriangularPoints()+getNumberOfEdges()*(degree-1)+triangleIndex*(degree-1)*(degree-2)/2+(*omi).second);
                 }
                 else
                 {
#ifndef NDEBUG
                    sout << "Unexpected error in [HighOrderTriangleSetTopologyContainer::getGlobalIndexOfControlPoint]" << sendl;
#endif
                    return 0; //Warning fix but maybe the author of this code would want to print a more meaningful error message for this "ei[0] = ei[1] = ei[2] = -1" case ?
                 }
             } else {
#ifndef NDEBUG
                 sout << "Error. [HighOrderTriangleSetTopologyContainer::getGlobalIndexOfControlPoint] Bezier Index "<< (sofa::defaulttype::Vec<3,int> )(id) <<" has not been recognized to be valid" << sendl;
#endif
                 return 0;
             }
         } else {
             std::map<ControlPointLocation,size_t>::const_iterator itgi;

             itgi=locationToGlobalIndexMap.find(ControlPointLocation(triangleIndex,id));
             if (itgi!=locationToGlobalIndexMap.end()) {
                 return(itgi->second);
             } else {
#ifndef NDEBUG
                 sout << "Error. [HighOrderTriangleSetTopologyContainer::getGlobalIndexOfControlPoint] Cannot find global index of control point with TRBI  "<< (sofa::defaulttype::Vec<3,int> )(id) <<" and triangle index " << triangleIndex <<sendl;
#endif
                 return 0;
             }

         }
 }


TriangleIndexVector HighOrderTriangleSetTopologyContainer::getTriangleIndexVector(const size_t localIndex) const
{
	
	if (localIndex<bezierIndexArray.size()) {
		return bezierIndexArray[localIndex];
	} else {
#ifndef NDEBUG
		sout << "Error. [HighOrderTriangleSetTopologyContainer::getTriangleIndexVector] Index "<< localIndex <<" is greater than the number "<< bezierIndexArray.size() <<" of control points." << sendl;
#endif
		TriangleIndexVector id;
		return (id);
	}
}
size_t HighOrderTriangleSetTopologyContainer::getLocalIndexFromTriangleIndexVector(const TriangleIndexVector id) const {
	OffsetMapConstIterator omi=localIndexMap.find(id);
	if (omi==localIndexMap.end())
	{
#ifndef NDEBUG
		sout << "Error. [HighOrderTriangleSetTopologyContainer::getLocalIndexFromTriangleIndexVector] Triangle Bezier Index "<< id  <<" is out of range." << sendl;
#endif
		return(0);
	} else {
		return ((*omi).second);
	}
}
size_t HighOrderTriangleSetTopologyContainer::getLexicographicIndex(size_t i) const {
	return(lexicographicIndexArray[i]);
}
size_t HighOrderTriangleSetTopologyContainer::getHierarchicalIndex(size_t i) const {
	return(hierarchicalIndexArray[i]);
}

size_t HighOrderTriangleSetTopologyContainer::getLexicographicIndex(const TriangleIndexVector tbi1) const {
	HighOrderDegreeType degree=d_degree.getValue();
	return((size_t)( (degree+1)*(degree+2)/2-(degree+1-tbi1[0])*(degree+2-tbi1[0])/2+tbi1[1]));
}
sofa::helper::vector<TriangleIndexVector> HighOrderTriangleSetTopologyContainer::getTriangleIndexArray() const
{
	return (bezierIndexArray);
}
sofa::helper::vector<TriangleIndexVector> HighOrderTriangleSetTopologyContainer::getTriangleIndexArrayOfGivenDegree(const HighOrderDegreeType deg) const
{
	// vertex index
	size_t i,j;
	sofa::helper::vector<TriangleIndexVector> tbiArray;
	for (i=0;i<3;++i) {
		TriangleIndexVector bti(0,0,0);
		bti[i]=deg;
		tbiArray.push_back(bti);
	}
	// edge index
	if (deg>1) {
		for (i=0;i<3;++i) {
			for (j=1;j<deg;++j) {
				TriangleIndexVector bti(0,0,0);
				bti[(i+1)%3]=(size_t)(deg-j);
				bti[(i+2)%3]=j;
				tbiArray.push_back(bti);
			}
		}
	}

	// Triangle index
	if (deg>2) {
		size_t ind=0;
        for (i=1;i<(HighOrderDegreeType)(deg-1);++i) {
			for (j=1;j<(deg-i);++j,++ind) {
				TriangleIndexVector bti(0,0,0);
				bti[0]=i;bti[1]=j;
				bti[2]=deg-i-j;
				tbiArray.push_back(bti);
			}
		}
	}

	return(tbiArray);
}
sofa::helper::vector<HighOrderTriangleSetTopologyContainer::LocalTriangleIndex> HighOrderTriangleSetTopologyContainer::getMapOfTriangleIndexArrayFromInferiorDegree() const
{
	HighOrderDegreeType degree=d_degree.getValue();
	sofa::helper::vector<TriangleIndexVector> tbiDerivArray=getTriangleIndexArrayOfGivenDegree(degree-1);
	sofa::helper::vector<TriangleIndexVector> tbiLinearArray=getTriangleIndexArrayOfGivenDegree(1);
	TriangleIndexVector tbi;
	sofa::helper::vector<LocalTriangleIndex> correspondanceArray;
//	correspondanceArray.resize(tbiDerivArray.size());
	size_t i,j;
	for (i=0;i<tbiDerivArray.size();++i) {
		LocalTriangleIndex correspondance;
		for (j=0;j<3;++j) {
			tbi=tbiDerivArray[i]+tbiLinearArray[j];
			correspondance[j]=getLocalIndexFromTriangleIndexVector(tbi);
		}
		correspondanceArray.push_back(correspondance);
	}
	return(correspondanceArray);
}
const HighOrderTriangleSetTopologyContainer::VecPointID &HighOrderTriangleSetTopologyContainer::getGlobalIndexArrayOfControlPoints(const TriangleID triangleIndex) const {
	return triangleDOFArray[triangleIndex];
}
void HighOrderTriangleSetTopologyContainer::getGlobalIndexArrayOfControlPointsInTriangle(const TriangleID triangleIndex, VecPointID& indexArray)
{
	
	indexArray.clear();
	if (locationToGlobalIndexMap.empty()) {
		Triangle tr=getTriangle(triangleIndex);
		// vertex index
		size_t i,j;
		for (i=0;i<3;++i)
			indexArray.push_back(tr[i]);

		size_t offset;
		// edge index
		HighOrderDegreeType degree=d_degree.getValue();
		if (degree>1) {
			EdgesInTriangle eit=getEdgesInTriangle(triangleIndex);
			for (i=0;i<3;++i) {
				Edge e=getEdge(eit[i]);
				offset=getNumberOfTriangularPoints()+eit[i]*(degree-1);
				// check the order of the edge to be consistent with the Triangle
				if (e[0]==tr[(i+1)%3] ) {
					for (j=0;j<(size_t)(degree-1);++j) {
						indexArray.push_back(offset+j);
					}
				} else {
					int jj;
					for (jj=degree-2;jj>=0;--jj) {
						indexArray.push_back(offset+jj);
					}
				}
			}
		}



		// Triangle index
		if (degree>2) {
			size_t pointsPerTriangle=(degree-1)*(degree-2)/2;
			offset=getNumberOfTriangularPoints()+getNumberOfEdges()*(degree-1)+triangleIndex*pointsPerTriangle;
			size_t rank=0;
			for (i=0;i<(size_t)(degree-2);++i) {
				for (j=0;j<(degree-i-2);++j) {
					indexArray.push_back(offset+rank);
					rank++;
				}
			}
		}
	} else {
		size_t i;
		std::map<ControlPointLocation,size_t>::const_iterator itgi;
		for (i=0;i<bezierIndexArray.size();++i) {
			itgi=locationToGlobalIndexMap.find(ControlPointLocation(triangleIndex,bezierIndexArray[i]));
			if (itgi!=locationToGlobalIndexMap.end()) {
				indexArray.push_back(itgi->second);
			} else {
#ifndef NDEBUG
				sout << "Error. [HighOrderTriangleSetTopologyContainer::getGlobalIndexArrayOfControlPointsInTriangle] Cannot find global index of control point with TRBI  "<< (sofa::defaulttype::Vec<3,int> )(bezierIndexArray[i]) <<" and triangle index " << triangleIndex <<sendl;
#endif
			}

		}

	}
	}
void HighOrderTriangleSetTopologyContainer::getLocationFromGlobalIndex(const size_t globalIndex, HighOrderTrianglePointLocation &location, 
	size_t &elementIndex, size_t &elementOffset)
{
	size_t gi=globalIndex;
	if (globalIndexToLocationMap.empty()) {
		if (gi<getNumberOfTriangularPoints()) {
			location=POINT;
			elementIndex=gi;
			elementOffset=0;
		} else {
			gi-=getNumberOfTriangularPoints();
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
#ifndef NDEBUG
					sout << "Error. [HighOrderTriangleSetTopologyContainer::getLocationFromGlobalIndex] Global Index "<< globalIndex <<" exceeds the number of Bezier Points" << sendl;
#endif
				}
			}

		}
	} else {
		std::multimap<size_t,ControlPointLocation>::const_iterator itcpl;
		itcpl=globalIndexToLocationMap.find(gi); 
		if (itcpl!=globalIndexToLocationMap.end()) {
			// get the local index and triangle index of that control point
			size_t offset=getLocalIndexFromTriangleIndexVector(itcpl->second.second);
			// if its local index is less than 3 then it is a triangle vertex
			if (offset<3) {
				location=POINT;
				elementIndex=getTriangle(itcpl->second.first)[offset];
				elementOffset=0;
			} else {
				offset -= 3;
				HighOrderDegreeType degree=d_degree.getValue();
                if ((HighOrderDegreeType)offset<3*(degree-1)){
					location=EDGE;
					// get the id of the edge on which it lies
					elementIndex=getEdgesInTriangle(itcpl->second.first)[offset/(degree-1)];
					elementOffset=offset%(degree-1);
				} else {
					offset -= 3*(degree-1);
					location=TRIANGLE;
					elementIndex=itcpl->second.first;
					elementOffset=offset;
				}
			}
		} else {
			location=NONE;
			elementIndex=0;
			elementOffset=0;
		}
	}
}
sofa::helper::vector<HighOrderTriangleSetTopologyContainer::LocalTriangleIndex> HighOrderTriangleSetTopologyContainer::getLocalIndexSubtriangleArray() const {
	sofa::helper::vector<LocalTriangleIndex> subtriangleArray;
	HighOrderDegreeType degree=d_degree.getValue();
	TriangleIndexVector tbi1,tbi2,tbi3;
	LocalTriangleIndex lti;
	for (size_t i=1;i<=degree;++i) {
		for (size_t j=0;j<(degree-i+1);++j) {
			tbi1=TriangleIndexVector(i,j,degree-i-j);
			tbi2=TriangleIndexVector(i-1,j+1,degree-i-j);
			tbi3=TriangleIndexVector(i-1,j,degree-i-j+1);
			
			lti[0]=getLocalIndexFromTriangleIndexVector(tbi1);
			lti[1]=getLocalIndexFromTriangleIndexVector(tbi2);
			lti[2]=getLocalIndexFromTriangleIndexVector(tbi3);
			subtriangleArray.push_back(lti);
			if ((i+j)<degree) {
				tbi3=TriangleIndexVector(i,j+1,degree-i-j-1);
				lti[2]=lti[1];	
				lti[1]=getLocalIndexFromTriangleIndexVector(tbi3);
				subtriangleArray.push_back(lti);
			}
		}
	}
    assert(subtriangleArray.size()==(size_t)((degree+1)*(degree+1)));
	return(subtriangleArray);

}
sofa::helper::vector<HighOrderTriangleSetTopologyContainer::LocalTriangleIndex> HighOrderTriangleSetTopologyContainer::getLocalIndexSubtriangleArrayOfGivenDegree(const HighOrderDegreeType deg)  const {

	sofa::helper::vector<TriangleIndexVector> tbia=getTriangleIndexArrayOfGivenDegree(deg);
	// create a local map for indexing
	std::map<TriangleIndexVector,size_t> tmpLocalIndexMap;
	size_t i;
	for (i=0;i<tbia.size();++i)
		tmpLocalIndexMap.insert(OffsetMapType(tbia[i],i));
	// now create the array of subtriangles
	sofa::helper::vector<LocalTriangleIndex> subtriangleArray;
	
	TriangleIndexVector tbi[3];
	size_t k;
	LocalTriangleIndex lti;
	std::map<TriangleIndexVector,size_t>::iterator omi;
	for ( i=1;i<=deg;++i) {
		for (size_t j=0;j<(deg-i+1);++j) {
			tbi[0]=TriangleIndexVector(i,j,deg-i-j);
			tbi[1]=TriangleIndexVector(i-1,j+1,deg-i-j);
			tbi[2]=TriangleIndexVector(i-1,j,deg-i-j+1);
			for (k=0;k<3;++k) {
				omi=tmpLocalIndexMap.find(tbi[k]);
				if (omi==tmpLocalIndexMap.end())
				{
#ifndef NDEBUG
					sout << "Error. [HighOrderTriangleSetTopologyContainer::getLocalIndexSubtriangleArrayOfGivenDegree(const HighOrderDegreeType deg) ] Triangle Bezier Index "<< tbi[k]  <<" is out of range." << sendl;
#endif
				} else {
					lti[k]= (*omi).second;
				}
			}

			subtriangleArray.push_back(lti);
			if ((i+j)<deg) {
				tbi[2]=TriangleIndexVector(i,j+1,deg-i-j-1);
				lti[2]=lti[1];
				omi=tmpLocalIndexMap.find(tbi[2]);
				if (omi==tmpLocalIndexMap.end())
				{
#ifndef NDEBUG
					sout << "Error. [HighOrderTriangleSetTopologyContainer::getLocalIndexSubtriangleArrayOfGivenDegree(const HighOrderDegreeType deg) ] Triangle Bezier Index "<< tbi[2]  <<" is out of range." << sendl;
#endif
				} else {
					lti[1]= (*omi).second;
				}

				subtriangleArray.push_back(lti);
			}
		}
	}
	assert(subtriangleArray.size()==((deg)*(deg)));
	return(subtriangleArray);
}
void HighOrderTriangleSetTopologyContainer::getEdgeIndexVectorFromEdgeOffset(size_t offset, EdgeIndexVector &ebi){
	assert(offset<d_degree.getValue());
	ebi[0]=offset+1;
	ebi[1]=d_degree.getValue()-offset-1;
}
void HighOrderTriangleSetTopologyContainer::getTriangleIndexVectorFromTriangleOffset(size_t offset, TriangleIndexVector &tbi){
	assert(offset<(d_degree.getValue()-1)*(d_degree.getValue()-2)/2);
	tbi=offsetToTriangleIndexVectorArray[offset];
}
bool HighOrderTriangleSetTopologyContainer::checkHighOrderTriangleTopology()
{
	#ifndef NDEBUG
	size_t nTrians,elem;
	HighOrderDegreeType degree=d_degree.getValue();
	// check the total number of vertices.
    assert(getNbPoints()==(int)(getNumberOfTriangularPoints()+getNumberOfEdges()*(degree-1)+getNumberOfTriangles()*(degree-1)*(degree-2)/2));
	sofa::helper::vector<TriangleIndexVector> tbiArray=getTriangleIndexArray();
	VecPointID indexArray;
	HighOrderTrianglePointLocation location; 
    size_t elementIndex, elementOffset/*,localIndex*/;
	for (nTrians=0;nTrians<getNumberOfTriangles();++nTrians) {
		indexArray.clear();
		getGlobalIndexArrayOfControlPointsInTriangle(nTrians,indexArray);
		// check the number of control points per Triangle is correct
        assert(indexArray.size()==(size_t)(3+3*(degree-1)+(degree-1)*(degree-2)/2));
		for(elem=0;elem<indexArray.size();++elem) {
			size_t globalIndex=getGlobalIndexOfControlPoint(nTrians,tbiArray[elem]);
			// check that getGlobalIndexOfControlPoint and getGlobalIndexArrayOfControlPointsInTriangle give the same answer
			assert(globalIndex==indexArray[elem]);

            TriangleIndexVector tbi=getTriangleIndexVector(elem);

			assert(elem==getLocalIndexFromTriangleIndexVector(tbi));
			// check that getTriangleIndexVector is consistent with getTriangleIndexVectorArray
			assert(tbiArray[elem][0]==tbi[0]);
			assert(tbiArray[elem][1]==tbi[1]);
			assert(tbiArray[elem][2]==tbi[2]);
			// check that getLocationFromGlobalIndex is consistent with 
			getLocationFromGlobalIndex(globalIndex,location,elementIndex,elementOffset);
			if (elem<3) {
				assert(location==POINT);
				assert(elementIndex==getTriangle(nTrians)[elem]);
				assert(elementOffset==0);
			}
			else if (elem<(size_t)(3+3*(degree-1))){
				assert(location==EDGE);
				assert(elementIndex==getEdgesInTriangle(nTrians)[(elem-3)/(degree-1)]);
			}
		}

	}
	if (locationToGlobalIndexMap.size()>0) {
		// check consistency between both maps
		assert(locationToGlobalIndexMap.size()==globalIndexToLocationMap.size());
		std::map<ControlPointLocation,size_t>::iterator itcpl;
		std::map<size_t,ControlPointLocation>::iterator itgi;
		std::pair<std::map<size_t,ControlPointLocation>::iterator,std::map<size_t,ControlPointLocation>::iterator> itgir;
		for (itcpl=locationToGlobalIndexMap.begin();itcpl!=locationToGlobalIndexMap.end();++itcpl) {
			itgir=globalIndexToLocationMap.equal_range(itcpl->second);
			assert(itgir.first!=itgir.second);
			for (itgi=itgir.first;itgi->second!=itcpl->first && itgi!=itgir.second;++itgi);
			assert(itgi->second==itcpl->first);
		}
	}
	#endif
	return( true);
}

} // namespace topology

} // namespace component

} // namespace sofa
