#ifndef SOFA_HIGHORDERTOPOLOGY_HIGHORDERTRIANGLESETGEOMETRYALGORITHMS_INL
#define SOFA_HIGHORDERTOPOLOGY_HIGHORDERTRIANGLESETGEOMETRYALGORITHMS_INL

#include "HighOrderTriangleSetGeometryAlgorithms.h"
#include <sofa/core/visual/VisualParams.h>
#include <SofaBaseTopology/CommonAlgorithms.h>
#include <sofa/helper/rmath.h>
namespace sofa
{

namespace component
{

namespace topology
{

template< class DataTypes>
 HighOrderTriangleSetGeometryAlgorithms< DataTypes >::HighOrderTriangleSetGeometryAlgorithms() : 
TriangleSetGeometryAlgorithms<DataTypes>()
		,degree(0)
        ,drawControlPointsEdges (core::objectmodel::Base::initData(&drawControlPointsEdges, (bool) false, "drawControlPointsEdges", "Debug : draw Control point edges "))
        ,drawSmoothEdges (core::objectmodel::Base::initData(&drawSmoothEdges, (bool) false, "drawSmoothEdges", "Debug : draw Bezier curves as edges of the  Bezier triangle"))
 	   ,drawControlPoints (core::objectmodel::Base::initData(&drawControlPoints, (bool) false, "drawControlPoints", "Debug : draw Control points with a color depending on its status "))
	   ,d_referenceRadius (core::objectmodel::Base::initData(&d_referenceRadius, (Real) -1.0f, "referenceRadius", "Debug : radius of control points when drawing "))
    {
    }
template< class DataTypes>
 void HighOrderTriangleSetGeometryAlgorithms< DataTypes >::init()
{
	TriangleSetGeometryAlgorithms<DataTypes>::init();
	// recovers the pointer to a BezierTriangleSetTopologyContainer
	HighOrderTriangleSetTopologyContainer *btstc = NULL;
	this->getContext()->get(btstc, sofa::core::objectmodel::BaseContext::SearchUp);
	if (!btstc) {
		serr << " Could not find a HighOrderTriangleSetTopologyContainer object"<< sendl;
	} else {
		container=btstc;
		/// get the degree of the Bezier tetrahedron
		degree=container->getDegree();
		/// store the tetrahedron bezier index for each tetrahedron
		tbiArray=container->getTriangleIndexArray();
		
		/// fills the array of edges
		bezierTriangleEdgeSet.clear();
		TriangleIndexVector tbiNext,tbi;
		size_t i,j,k,m,n,index1,index2;
		for (i=0;i<=degree;++i) {
			for (j=0;j<=(degree-i);++j) {
				k=degree-i-j;
				tbi=TriangleIndexVector(i,j,k);
				index1=container->getLocalIndexFromTriangleIndexVector(tbi);
				for(m=0;m<3;++m) {
					if (tbi[m]<degree) {
						for (n=1;n<3;++n) {
							if (tbi[(m+n)%3]!=0) {
								tbiNext=tbi;
								tbiNext[m]=tbi[m]+1;
								tbiNext[(m+n)%3]=tbi[(m+n)%3]-1;
								index2=container->getLocalIndexFromTriangleIndexVector(tbiNext);
								Edge e((PointID)std::min(index1,index2),(PointID)std::max(index1,index2));
								// test if both control points are on an edge
								if (tbi[(m+3-n)%3]==0)
									bezierTriangleEdgeSet.insert(std::pair<Edge,bool>(e,true));
								else 
									bezierTriangleEdgeSet.insert(std::pair<Edge,bool>(e,false));
							}
						}
					}
				}
			}
		}

	}

	if (this->f_printLog.getValue()) {
		size_t nbAffine=0;
		const typename DataTypes::VecCoord& p =(this->object->read(core::ConstVecCoordId::position())->getValue());
		for (size_t i=0;i<(size_t)container->getNbTriangles();++i) {
			if (isBezierTriangleAffine(i,p))
				nbAffine++;
		}
		sout << "Number of affine triangles   :"<< nbAffine<< " / "<< (size_t)container->getNbTriangles() <<sendl;
	}
 }

template< class DataTypes>
 void HighOrderTriangleSetGeometryAlgorithms< DataTypes >::reinit()
{
}

template< class DataTypes>
typename DataTypes::Coord HighOrderTriangleSetGeometryAlgorithms< DataTypes >::computeNodalValue(const size_t triangleIndex,const Vec3 barycentricCoordinate, const typename DataTypes::VecCoord& p)
{
	Coord nodalValue;
	nodalValue.clear();

	TriangleIndexVector tbi;
	const VecPointID &indexArray=container->getGlobalIndexArrayOfControlPoints(triangleIndex);
	bool isRational=container->isRationalSpline(triangleIndex);
	if (isRational) {
		const HighOrderTriangleSetTopologyContainer::SeqWeights &wa=container->getWeightArray();
		Real weight=(Real)0.0f;
		Real shapeFunctionValue;
		for(size_t i=0; i<tbiArray.size(); ++i)
		{
			tbi=tbiArray[i];
			shapeFunctionValue=computeShapeFunction(tbi,barycentricCoordinate); 
			nodalValue+=wa[indexArray[i]]*p[indexArray[i]]*shapeFunctionValue;
			weight+=wa[indexArray[i]]*shapeFunctionValue;
		}
		nodalValue/=weight;
	} else {
		for(size_t i=0; i<tbiArray.size(); ++i)
		{
			tbi=tbiArray[i];
			nodalValue+=p[indexArray[i]]*computeShapeFunction(tbi,barycentricCoordinate); 
		}
	}
	return(nodalValue);
}
template< class DataTypes>
typename DataTypes::Coord HighOrderTriangleSetGeometryAlgorithms< DataTypes >::computeNodalValue(const size_t triangleIndex,const Vec3 barycentricCoordinate)
{
	const typename DataTypes::VecCoord& p =(this->object->read(core::ConstVecCoordId::position())->getValue());
	return(computeNodalValue(triangleIndex,barycentricCoordinate,p));
}


 template<class DataTypes>
 typename HighOrderTriangleSetGeometryAlgorithms<DataTypes>::Real 
	 HighOrderTriangleSetGeometryAlgorithms<DataTypes>::computeJacobian(const size_t triangleIndex, const Vec3 barycentricCoordinate)
 {
	 const typename DataTypes::VecCoord& p =(this->object->read(core::ConstVecCoordId::position())->getValue());
	 return(computeJacobian(triangleIndex,barycentricCoordinate,p));

 }

template<class DataTypes>
bool HighOrderTriangleSetGeometryAlgorithms<DataTypes>::isBezierTriangleAffine(const size_t triangleIndex,const VecCoord& p, Real tolerance) const{
	// get the global indices of all points
	const VecPointID &indexArray=container->getGlobalIndexArrayOfControlPoints(triangleIndex);
	bool affine=true;

	/// skip the first 4 control points corresponding to the 4 corners
	size_t index=0;
	Coord corner[3],pos,actualPos;
	// store the position of the 4 corners
	for (index=0;index<3;++index)
		corner[index]=p[indexArray[index]];
	do {
		// compute the position of the control point as if the Triangle was affine
		pos=corner[0]*tbiArray[index][0]+corner[1]*tbiArray[index][1]+corner[2]*tbiArray[index][2];
			
		pos/=degree;
		// measure the distance between the real position and the affine position
		actualPos=p[indexArray[index]];
		if ((actualPos-pos).norm2()>tolerance) {
			affine=false;
		}
		index++;
	} while ((affine) && (index<indexArray.size()));
	return (affine);
}

template<class DataTypes>
void HighOrderTriangleSetGeometryAlgorithms<DataTypes>::computeNodalValueDerivatives(const size_t triangleIndex, const Vec3 barycentricCoordinate,  const VecCoord& p, Deriv dpos[3])
{
	/// the 4 derivatives
	size_t j;
	const VecPointID &indexArray=container->getGlobalIndexArrayOfControlPoints(triangleIndex);
	// initialize dpos
	for (j=0;j<3;++j) 
		dpos[j]=Deriv();
	for(size_t i=0; i<tbiArray.size(); ++i)
	{
		Vec3 dval=computeShapeFunctionDerivatives(tbiArray[i], barycentricCoordinate);
		for (j=0;j<3;++j) {
			dpos[j]+=dval[j]*p[indexArray[i]];
		}
	}
}

template<class DataTypes>
typename HighOrderTriangleSetGeometryAlgorithms<DataTypes>::Real 
	HighOrderTriangleSetGeometryAlgorithms<DataTypes>::computeJacobian(const size_t triangleIndex, const Vec3 barycentricCoordinate, const typename DataTypes::VecCoord& p)
{
	/// the 3 derivatives
	Deriv dpos[3];

	computeNodalValueDerivatives(triangleIndex,barycentricCoordinate,p,dpos);
	return(areaProduct(dpos[0]-dpos[2],dpos[1]-dpos[2]));
}





template<class DataTypes>
void HighOrderTriangleSetGeometryAlgorithms<DataTypes>::draw(const core::visual::VisualParams* vparams)
{

	
	if ((degree>0) && (container) ){
		TriangleSetGeometryAlgorithms<DataTypes>::draw(vparams);	
		if (drawControlPoints.getValue())
		{
			size_t nbPoints=container->getNbPoints();
			size_t i,elementIndex,elementOffset;
			const typename DataTypes::VecCoord& pos =(this->object->read(core::ConstVecCoordId::position())->getValue());
			HighOrderTriangleSetTopologyContainer::HighOrderTrianglePointLocation location;

			if (container->getNbTriangles()>0) {
				// estimate the  mean radius of the spheres from the first Bezier triangle 
				const VecPointID &indexArray=container->getGlobalIndexArrayOfControlPoints(0);
				std::vector<Real> edgeLengthArray;
				// compute median of the edge distance between control points	
				std::set<std::pair<Edge,size_t> >::iterator ite=bezierTriangleEdgeSet.begin();
//				Real val=0;
				Coord pp;
				for (; ite!=bezierTriangleEdgeSet.end(); ite++)
				{
						pp = pos[indexArray[(*ite).first[0]]] -pos[indexArray[(*ite).first[1]]] ;
						edgeLengthArray.push_back(pp.norm());
				}
				std::nth_element(edgeLengthArray.begin(), edgeLengthArray.begin() + edgeLengthArray.size()/2, edgeLengthArray.end());
				Real radius=edgeLengthArray[edgeLengthArray.size()/2]/5;
				std::vector<sofa::defaulttype::Vector3> pointsVertices,pointsEdges,pointsTriangles;
				std::vector<float> radiusVertices,radiusEdges,radiusTriangles;
				sofa::defaulttype::Vector3 p;


				for (i=0;i<nbPoints;++i) {
					container->getLocationFromGlobalIndex(i,location,elementIndex,elementOffset);
					if (location==HighOrderTriangleSetTopologyContainer::NONE) {
					} else if (location==HighOrderTriangleSetTopologyContainer::POINT) {
						p=pos[i];
						pointsVertices.push_back(p);

						radiusVertices.push_back(radius*container->getWeight(i));

					} else if (location==HighOrderTriangleSetTopologyContainer::EDGE) {
						p=pos[i];
						pointsEdges.push_back(p);

						radiusEdges.push_back(radius*container->getWeight(i));

					} else {
						p=pos[i];
						pointsTriangles.push_back(p);

						radiusTriangles.push_back(radius*container->getWeight(i));

					}
				}
				vparams->drawTool()->setLightingEnabled(true); //Enable lightning
				vparams->drawTool()->drawSpheres(pointsVertices, radiusVertices,  defaulttype::Vec<4,float>(1.0f,0,0,1.0f));
				vparams->drawTool()->drawSpheres(pointsEdges, radiusEdges,  defaulttype::Vec<4,float>(0,1.0f,0,1.0f));
				vparams->drawTool()->drawSpheres(pointsTriangles, radiusTriangles,  defaulttype::Vec<4,float>(0,0,1.0f,1.0f));
				vparams->drawTool()->setLightingEnabled(false); //Disable lightning
			}
		}
		// Draw edges linking Bezier Triangle control points with a color code
		if (drawSmoothEdges.getValue())
		{
			
			const sofa::helper::vector<Triangle> &trianArray = this->m_topology->getTriangles();

			if (!trianArray.empty())
			{
				// estimate the  mean radius of the spheres from the first Bezier triangle 
				size_t nbPoints=container->getNbPoints();
				size_t i,elementIndex,elementOffset;
				const typename DataTypes::VecCoord& coords =(this->object->read(core::ConstVecCoordId::position())->getValue());
				HighOrderTriangleSetTopologyContainer::HighOrderTrianglePointLocation location;
				const VecPointID &indexArrayFirst=container->getGlobalIndexArrayOfControlPoints(0);
				std::vector<Real> edgeLengthArray;
				// compute median of the edge distance between control points	
				std::set<std::pair<Edge,size_t> >::iterator ite=bezierTriangleEdgeSet.begin();
//				Real val=0;
				Coord pp;
				for (; ite!=bezierTriangleEdgeSet.end(); ite++)
				{
					pp = coords[indexArrayFirst[(*ite).first[0]]] -coords[indexArrayFirst[(*ite).first[1]]] ;
					edgeLengthArray.push_back(pp.norm());
				}
				std::nth_element(edgeLengthArray.begin(), edgeLengthArray.begin() + edgeLengthArray.size()/2, edgeLengthArray.end());
				Real radius=edgeLengthArray[edgeLengthArray.size()/2]/5;
				std::vector<sofa::defaulttype::Vector3> pointsVertices;
				std::vector<float> radiusVertices;
				sofa::defaulttype::Vector3 p1;


				for (i=0;i<nbPoints;++i) {
					container->getLocationFromGlobalIndex(i,location,elementIndex,elementOffset);
					if (location==HighOrderTriangleSetTopologyContainer::POINT) {
						p1=coords[i];
						pointsVertices.push_back(p1);

						radiusVertices.push_back(radius*container->getWeight(i));

					} 
				}
				vparams->drawTool()->setLightingEnabled(true); //Enable lightning
				vparams->drawTool()->drawSpheres(pointsVertices, radiusVertices,  defaulttype::Vec<4,float>(1.0f,0,0,1.0f));
				vparams->drawTool()->setLightingEnabled(false); //Disable lightning

				#ifndef SOFA_NO_OPENGL
				glDisable(GL_LIGHTING);
				
				glColor3f(0.0f, 1.0f, 0.0f);
				glLineWidth(3.0);
				 glEnable(GL_DEPTH_TEST);
				glEnable(GL_POLYGON_OFFSET_LINE);
				glPolygonOffset(-1.0,100.0);
				
				
				// how many points is used to discretize the edge
				const size_t edgeTesselation=9;
				sofa::defaulttype::Vec3f p; //,p2;
				for ( i = 0; i<trianArray.size(); i++)
				{
					const VecPointID &indexArray=container->getGlobalIndexArrayOfControlPoints(i);
					sofa::helper::vector <sofa::defaulttype::Vec3f> trianCoord;
					// process each edge of the triangle
					for (size_t j = 0; j<3; j++) {
						Vec3 baryCoord;
						baryCoord[j]=0;
						glBegin(GL_LINE_STRIP);
						for (size_t k=0;k<=edgeTesselation;++k) {
							baryCoord[(j+1)%3]=(Real)k/(Real)edgeTesselation;
							baryCoord[(j+2)%3]=(Real)(edgeTesselation-k)/(Real)edgeTesselation;
							p=DataTypes::getCPos(computeNodalValue(i,baryCoord));
							glVertex3f(p[0],p[1],p[2]);
						}
						glEnd();
					}
				}
				glDisable(GL_POLYGON_OFFSET_LINE);
				
			}
			#endif
		}
		
		// Draw edges linking Bezier Triangle control points with a color code
		if (drawControlPointsEdges.getValue())
		{
			#ifndef SOFA_NO_OPENGL
			const sofa::helper::vector<Triangle> &trianArray = this->m_topology->getTriangles();

			if (!trianArray.empty())
			{
				glDisable(GL_LIGHTING);
                const sofa::defaulttype::Vec4f& color =  this->_drawColor.getValue();
				glColor3f(color[0], color[1], color[2]);
				glBegin(GL_LINES);
				const VecCoord& coords =(this->object->read(core::ConstVecCoordId::position())->getValue());
				VecPointID indexArray;
				Vec3 baryCoord;
				sofa::defaulttype::Vec3f p; //,p2;
				for (unsigned int i = 0; i<trianArray.size(); i++)
				{
					const VecPointID &indexArray=container->getGlobalIndexArrayOfControlPoints(i);
					sofa::helper::vector <sofa::defaulttype::Vec3f> trianCoord;

					for (unsigned int j = 0; j<indexArray.size(); j++)
					{
						p = DataTypes::getCPos(coords[indexArray[j]]);
						trianCoord.push_back(p);
					}
					
					std::set<std::pair<Edge,size_t> >::iterator ite=bezierTriangleEdgeSet.begin();
					for (; ite!=bezierTriangleEdgeSet.end(); ite++)
					{
						if ((*ite).second) {
							glColor3f(0.0f, 1.0f, 0.0f);
						} else {
							glColor3f(0.0f, 0.0f, 1.0f );
						}
						glVertex3f(trianCoord[(*ite).first[0]][0], trianCoord[(*ite).first[0]][1], trianCoord[(*ite).first[0]][2]);
						glVertex3f(trianCoord[(*ite).first[1]][0], trianCoord[(*ite).first[1]][1], trianCoord[(*ite).first[1]][2]);
					}
				}
				glEnd();
			}
			#endif
		}
	}
	

}



} // namespace topology

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENTS_HIGHORDERTETEAHEDRONSETGEOMETRYALGORITHMS_INL
