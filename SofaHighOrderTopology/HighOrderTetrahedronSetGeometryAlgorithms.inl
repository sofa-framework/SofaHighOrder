#ifndef SOFA_HIGHORDERTOPOLOGY_HIGHORDERTETRAHEDRONSETGEOMETRYALGORITHMS_INL
#define SOFA_HIGHORDERTOPOLOGY_HIGHORDERTETRAHEDRONSETGEOMETRYALGORITHMS_INL

#include "HighOrderTetrahedronSetGeometryAlgorithms.h"
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
 HighOrderTetrahedronSetGeometryAlgorithms< DataTypes >::HighOrderTetrahedronSetGeometryAlgorithms() : 
TetrahedronSetGeometryAlgorithms<DataTypes>()
		,degree(0)
        ,drawControlPointsEdges (core::objectmodel::Base::initData(&drawControlPointsEdges, (bool) false, "drawControlPointsEdges", "Debug : draw Control point edges "))
        ,drawSmoothEdges (core::objectmodel::Base::initData(&drawSmoothEdges, (bool) false, "drawSmoothEdges", "Debug : draw Bezier curves as edges of the  Bezier triangle"))
 	   ,drawControlPoints (core::objectmodel::Base::initData(&drawControlPoints, (bool) false, "drawControlPoints", "Debug : draw Control points with a color depending on its status "))
	   ,d_referenceRadius (core::objectmodel::Base::initData(&d_referenceRadius, (Real) -1.0f, "referenceRadius", "Debug : radius of control points when drawing "))
    {
    }
template< class DataTypes>
 void HighOrderTetrahedronSetGeometryAlgorithms< DataTypes >::init()
{
	TetrahedronSetGeometryAlgorithms<DataTypes>::init();
	// recovers the pointer to a BezierTetrahedronSetTopologyContainer
	HighOrderTetrahedronSetTopologyContainer *btstc = NULL;
	this->getContext()->get(btstc, sofa::core::objectmodel::BaseContext::SearchUp);
	if (!btstc) {
		serr << " Could not find a HighOrderTetrahedronSetTopologyContainer object"<< sendl;
	} else {
		container=btstc;
		/// get the degree of the Bezier tetrahedron
		degree=container->getDegree();
		/// store the tetrahedron bezier index for each tetrahedron
		tbiArray=container->getTetrahedronIndexArray();
		
		/// fills the array of edges
		bezierTetrahedronEdgeSet.clear();
		TetrahedronIndexVector tbiNext,tbi;
		/// insert coefficient for the inferior degree
		HighOrderDegreeType i,j,k,l,m,n,index1,index2;

		for (i=0;i<=degree;++i) {
			for (j=0;j<=(degree-i);++j) {
				for (k=0;k<=(degree-j-i);++k) {
					l=degree-i-j-k;
					tbi=TetrahedronIndexVector(i,j,k,l);
					index1=container->getLocalIndexFromTetrahedronIndexVector(tbi);
					for(m=0;m<4;++m) {
						if (tbi[m]<degree) {
							for (n=1;n<4;++n) {
								if (tbi[(m+n)%4]!=0) {
									tbiNext=tbi;
									tbiNext[m]=tbi[m]+1;
									tbiNext[(m+n)%4]=tbi[(m+n)%4]-1;
									index2=container->getLocalIndexFromTetrahedronIndexVector(tbiNext);
									Edge e((PointID)std::min(index1,index2),(PointID)std::max(index1,index2));
									// test if both control points are on an edge or an
									if (tbi[(m+1+(n%3))%4]==0) {
										if (tbi[(m+1+((n+1)%3))%4]==0) {
											// edge connects points along an edge
											bezierTetrahedronEdgeSet.insert(std::pair<Edge,size_t>(e,(size_t)2));
										} else 
											// edge connects points along a triangle
											bezierTetrahedronEdgeSet.insert(std::pair<Edge,size_t>(e,(size_t)1));
									} else if  (tbi[(m+1+((n+1)%3))%4]==0) {
										bezierTetrahedronEdgeSet.insert(std::pair<Edge,size_t>(e,(size_t)1));
									} else
										bezierTetrahedronEdgeSet.insert(std::pair<Edge,size_t>(e,(size_t)0));
								}
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
		for (size_t i=0;i<(size_t)container->getNbTetrahedra();++i) {
			if (isBezierTetrahedronAffine(i,p))
				nbAffine++;
		}
		sout << "Number of affine tetrahedra   :"<< nbAffine<< " / "<< (size_t)container->getNbTetrahedra() <<sendl;
	}
}

template< class DataTypes>
 void HighOrderTetrahedronSetGeometryAlgorithms< DataTypes >::reinit()
{
}

template< class DataTypes>
typename DataTypes::Coord HighOrderTetrahedronSetGeometryAlgorithms< DataTypes >::computeNodalValue(const size_t tetrahedronIndex,const Vec4 barycentricCoordinate, const typename DataTypes::VecCoord& p)
{
	Coord nodalValue;
	nodalValue.clear();

	TetrahedronIndexVector tbi;
	const VecPointID &indexArray=container->getGlobalIndexArrayOfControlPoints(tetrahedronIndex);
	bool isRational=container->isRationalSpline(tetrahedronIndex);
	if (isRational) {
		const HighOrderTetrahedronSetTopologyContainer::SeqWeights &wa=container->getWeightArray();
		Real weight=(Real)0.0f;
		Real bernsteinPolynonial;
		for(size_t i=0; i<tbiArray.size(); ++i)
		{
			tbi=tbiArray[i];
			bernsteinPolynonial=computeShapeFunction(tbi,barycentricCoordinate); 
			nodalValue+=wa[indexArray[i]]*p[indexArray[i]]*bernsteinPolynonial;
			weight+=wa[indexArray[i]]*bernsteinPolynonial;
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
typename DataTypes::Coord HighOrderTetrahedronSetGeometryAlgorithms< DataTypes >::computeNodalValue(const size_t tetrahedronIndex,const Vec4 barycentricCoordinate)
{
	const typename DataTypes::VecCoord& p =(this->object->read(core::ConstVecCoordId::position())->getValue());
	return(computeNodalValue(tetrahedronIndex,barycentricCoordinate,p));
}


 template<class DataTypes>
 typename HighOrderTetrahedronSetGeometryAlgorithms<DataTypes>::Real 
	 HighOrderTetrahedronSetGeometryAlgorithms<DataTypes>::computeJacobian(const size_t tetrahedronIndex, const Vec4 barycentricCoordinate)
 {
	 const typename DataTypes::VecCoord& p =(this->object->read(core::ConstVecCoordId::position())->getValue());
	 return(computeJacobian(tetrahedronIndex,barycentricCoordinate,p));

 }

template<class DataTypes>
bool HighOrderTetrahedronSetGeometryAlgorithms<DataTypes>::isBezierTetrahedronAffine(const size_t tetrahedronIndex,const VecCoord& p, Real tolerance) const{
	// get the global indices of all points
	
	const VecPointID &indexArray=container->getGlobalIndexArrayOfControlPoints(tetrahedronIndex);
	bool affine=true;

	/// skip the first 4 control points corresponding to the 4 corners
	size_t index=0;
	Coord corner[4],pos,actualPos;
	// store the position of the 4 corners
	for (index=0;index<4;++index)
		corner[index]=p[indexArray[index]];
	do {
		// compute the position of the control point as if the tetrahedron was affine
		pos=corner[0]*tbiArray[index][0]+corner[1]*tbiArray[index][1]+corner[2]*tbiArray[index][2]+
			corner[3]*tbiArray[index][3];
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
void HighOrderTetrahedronSetGeometryAlgorithms<DataTypes>::computeNodalValueDerivatives(const size_t tetrahedronIndex, const Vec4 barycentricCoordinate,  const VecCoord& p, Deriv dpos[4])
{
	/// the 4 derivatives
	size_t j;
	const VecPointID &indexArray=container->getGlobalIndexArrayOfControlPoints(tetrahedronIndex);
	// initialize dpos
	for (j=0;j<4;++j) 
		dpos[j]=Deriv();
	for(size_t i=0; i<tbiArray.size(); ++i)
	{
		Vec4 dval=computeShapeFunctionDerivatives(tbiArray[i], barycentricCoordinate);
		for (j=0;j<4;++j) {
			dpos[j]+=dval[j]*p[indexArray[i]];
		}
	}
}

template<class DataTypes>
typename HighOrderTetrahedronSetGeometryAlgorithms<DataTypes>::Real 
	HighOrderTetrahedronSetGeometryAlgorithms<DataTypes>::computeJacobian(const size_t tetrahedronIndex, const Vec4 barycentricCoordinate, const typename DataTypes::VecCoord& p)
{
	/// the 4 derivatives
	Deriv dpos[4];

	computeNodalValueDerivatives(tetrahedronIndex,barycentricCoordinate,p,dpos);
	return(tripleProduct(dpos[0]-dpos[3],dpos[1]-dpos[3],dpos[2]-dpos[3]));
}


template<class DataTypes>
void HighOrderTetrahedronSetGeometryAlgorithms<DataTypes>::draw(const core::visual::VisualParams* vparams)
{

	if ((degree>0) && (container) ){
		TetrahedronSetGeometryAlgorithms<DataTypes>::draw(vparams);	
		// Draw Tetra
		// reference radius
		if (d_referenceRadius.getValue()<0.0) {
			// estimate the  mean radius of the spheres from the first Bezier triangle 
			
//			size_t nbPoints=container->getNbPoints();
//			size_t i,elementIndex,elementOffset;
			const typename DataTypes::VecCoord& coords =(this->object->read(core::ConstVecCoordId::position())->getValue());
//			BezierTetrahedronSetTopologyContainer::BezierTetrahedronPointLocation location;
			const VecPointID &indexArray=container->getGlobalIndexArrayOfControlPoints(0);
			std::vector<Real> edgeLengthArray;
			// compute median of the edge distance between control points	
			std::set<std::pair<Edge,size_t> >::iterator ite=bezierTetrahedronEdgeSet.begin();
			//Real val=0;
			Coord pp;
			for (; ite!=bezierTetrahedronEdgeSet.end(); ite++)
			{
				pp = coords[indexArray[(*ite).first[0]]] -coords[indexArray[(*ite).first[1]]] ;
				edgeLengthArray.push_back(pp.norm());
			}
			std::nth_element(edgeLengthArray.begin(), edgeLengthArray.begin() + edgeLengthArray.size()/2, edgeLengthArray.end());
			Real radius=edgeLengthArray[edgeLengthArray.size()/2]/5;
			d_referenceRadius.setValue(radius);
		}
		
		if (drawControlPoints.getValue())
		{
			size_t nbPoints=container->getNbPoints();
			size_t i,elementIndex,elementOffset;
			const typename DataTypes::VecCoord& pos =(this->object->read(core::ConstVecCoordId::position())->getValue());
			HighOrderTetrahedronSetTopologyContainer::HighOrderTetrahedronPointLocation location;

			if (container->getNbTriangles()>0) {
				// estimate the  mean radius of the spheres from the first Bezier triangle 
				VecPointID indexArray;
				float radius=	d_referenceRadius.getValue();
				std::vector<sofa::defaulttype::Vector3> pointsVertices,pointsEdges,pointsTriangles,pointsTetrahedra;
				std::vector<float> radiusVertices,radiusEdges,radiusTriangles,radiusTetrahedra;
				sofa::defaulttype::Vector3 p;


				for (i=0;i<nbPoints;++i) {
					container->getLocationFromGlobalIndex(i,location,elementIndex,elementOffset);
					if (location==HighOrderTetrahedronSetTopologyContainer::POINT) {
						p=pos[i];
						pointsVertices.push_back(p);

						radiusVertices.push_back(radius*container->getWeight(i));

					} else if (location==HighOrderTetrahedronSetTopologyContainer::EDGE) {
						p=pos[i];
						pointsEdges.push_back(p);

						radiusEdges.push_back(radius*container->getWeight(i));

					} else if (location==HighOrderTetrahedronSetTopologyContainer::TRIANGLE) {
						p=pos[i];
						pointsTriangles.push_back(p);

						radiusTriangles.push_back(radius*container->getWeight(i));

					} else {
						p=pos[i];
						pointsTetrahedra.push_back(p);

						radiusTetrahedra.push_back(radius*container->getWeight(i));
					}
				}
				vparams->drawTool()->setLightingEnabled(true); //Enable lightning
				vparams->drawTool()->drawSpheres(pointsVertices, radiusVertices,  defaulttype::Vec<4,float>(1.0f,0,0,1.0f));
				vparams->drawTool()->drawSpheres(pointsEdges, radiusEdges,  defaulttype::Vec<4,float>(0,1.0f,0,1.0f));
				vparams->drawTool()->drawSpheres(pointsTriangles, radiusTriangles,  defaulttype::Vec<4,float>(0,0,1.0f,1.0f));
				vparams->drawTool()->drawSpheres(pointsTetrahedra, radiusTetrahedra,  defaulttype::Vec<4,float>(1,0,1.0f,1.0f));
				vparams->drawTool()->setLightingEnabled(false); //Disable lightning
			}
		}
		// Draw edges linking Bezier tetrahedra control points with a color code
		if (drawSmoothEdges.getValue())
		{
			
			const sofa::helper::vector<Tetrahedron> &tetraArray = this->m_topology->getTetrahedra();

			if (!tetraArray.empty())
			{
				float radius=	d_referenceRadius.getValue();
				const typename DataTypes::VecCoord& coords =(this->object->read(core::ConstVecCoordId::position())->getValue());
				std::vector<sofa::defaulttype::Vector3> pointsVertices;
				std::vector<float> radiusVertices;
				sofa::defaulttype::Vector3 p1;
				size_t elementIndex,i,elementOffset;
				HighOrderTetrahedronSetTopologyContainer::HighOrderTetrahedronPointLocation location;
				size_t nbPoints=container->getNbPoints();

				for (i=0;i<nbPoints;++i) {
					container->getLocationFromGlobalIndex(i,location,elementIndex,elementOffset);
					if (location==HighOrderTetrahedronSetTopologyContainer::POINT) {
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
				
				const unsigned int edgesInTetrahedronArray[6][2] = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}};
				const unsigned int oppositeEdgesInTetrahedronArray[6][2] = {{2,3}, {1,3}, {1,2}, {0,3}, {0,2}, {0,1}};
				// how many points is used to discretize the edge
				const size_t edgeTesselation=9;
				sofa::defaulttype::Vec3f p; //,p2;
				for ( i = 0; i<tetraArray.size(); i++)
				{
					
//					const VecPointID &indexArray=container->getGlobalIndexArrayOfBezierPoints(i);
//					sofa::helper::vector <sofa::defaulttype::Vec3f> trianCoord;
					// process each edge of the tetrahedron
					for (size_t j = 0; j<6; j++) {
						Vec4 baryCoord;
						baryCoord[oppositeEdgesInTetrahedronArray[j][0]]=0;
						baryCoord[oppositeEdgesInTetrahedronArray[j][1]]=0;
						glBegin(GL_LINE_STRIP);
						for (size_t k=0;k<=edgeTesselation;++k) {
							baryCoord[edgesInTetrahedronArray[j][0]]=(Real)k/(Real)edgeTesselation;
							baryCoord[edgesInTetrahedronArray[j][1]]=(Real)(edgeTesselation-k)/(Real)edgeTesselation;
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

		if (drawControlPointsEdges.getValue())
		{
			const sofa::helper::vector<Tetrahedron> &tetraArray = this->m_topology->getTetrahedra();
				#ifndef SOFA_NO_OPENGL
			if (!tetraArray.empty())
			{
				glDisable(GL_LIGHTING);
				const sofa::defaulttype::Vec4f& color =  this->_drawColor.getValue();
				glColor3f(color[0], color[1], color[2]);
				glBegin(GL_LINES);
				const VecCoord& coords =(this->object->read(core::ConstVecCoordId::position())->getValue());
				
				Vec4 baryCoord;
				sofa::defaulttype::Vec3f p; //,p2;
				for (unsigned int i = 0; i<tetraArray.size(); i++)
				{
					
					const VecPointID &indexArray=container->getGlobalIndexArrayOfControlPoints(i);
					sofa::helper::vector <sofa::defaulttype::Vec3f> tetraCoord;

					for (unsigned int j = 0; j<indexArray.size(); j++)
					{
						p = DataTypes::getCPos(coords[indexArray[j]]);
						tetraCoord.push_back(p);
					}

			
					std::set<std::pair<Edge,size_t> >::iterator ite=bezierTetrahedronEdgeSet.begin();
					for (; ite!=bezierTetrahedronEdgeSet.end(); ite++)
					{
						if ((*ite).second==2) {
							glColor3f(0.0f, 1.0f, 0.0f);
						} else 	if ((*ite).second==1)  {
							glColor3f(0.0f, 0.0f, 1.0f );
						} else {
							glColor3f(1.0f, 0.0f, 1.0f );
						}
						glVertex3f(tetraCoord[(*ite).first[0]][0], tetraCoord[(*ite).first[0]][1], tetraCoord[(*ite).first[0]][2]);
						glVertex3f(tetraCoord[(*ite).first[1]][0], tetraCoord[(*ite).first[1]][1], tetraCoord[(*ite).first[1]][2]);

						
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
