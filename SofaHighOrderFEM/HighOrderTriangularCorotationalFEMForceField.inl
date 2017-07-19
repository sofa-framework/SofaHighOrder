
#ifndef SOFA_COMPONENT_FORCEFIELD_HIGHORDERTRIANGULARCOROTATIONALFEMFORCEFIELD_INL
#define SOFA_COMPONENT_FORCEFIELD_HIGHORDERTRIANGULARCOROTATIONALFEMFORCEFIELD_INL

#include "HighOrderTriangularCorotationalFEMForceField.h"
#include <sofa/core/visual/VisualParams.h>
#include <fstream> // for reading the file
#include <iostream> //for debugging
#include <sofa/helper/gl/template.h>
#include <SofaBaseTopology/TopologyData.inl>
#include <HighOrderTriangleSetGeometryAlgorithms.h>
#include <BezierTriangleSetGeometryAlgorithms.h>
#include <sofa/core/behavior/ForceField.inl>
#include <SofaBaseTopology/CommonAlgorithms.h>
#include <sofa/helper/decompose.h>
#include <boost/make_shared.hpp>

namespace sofa
{

namespace component
{

namespace forcefield
{

using namespace sofa::defaulttype;
using namespace	sofa::component::topology;
using namespace core::topology;

using core::topology::BaseMeshTopology;

typedef BaseMeshTopology::Triangle				Triangle;
typedef BaseMeshTopology::EdgesInTriangle		EdgesInTriangle;

typedef Triangle			        Triangle;
typedef EdgesInTriangle		EdgesInTriangle;

const unsigned int edgesInTriangleArray[3][2] = {{0,1}, {0,2}, {1,2}};




template< class DataTypes>
void HighOrderTriangularCorotationalFEMForceField<DataTypes>::weightArrayPointer::allocate(size_t degree) {
	if (degree==2) {
		weightArrayQuadratic[0]=boost::make_shared<Mat15x3>();
		weightArrayQuadratic[1]=boost::make_shared<Mat15x3>();

	} else 	if (degree==3) {
		weightArrayCubic[0]=boost::make_shared<Mat45x3>();
		weightArrayCubic[1]=boost::make_shared<Mat45x3>();
	} else 	if (degree==4) {
		weightArrayQuartic[0]=boost::make_shared<Mat105x3>();
		weightArrayQuartic[1]=boost::make_shared<Mat105x3>();
	} else 	if (degree==5) {
		weightArrayQuintic[0]=boost::make_shared<Mat210x3>();
		weightArrayQuintic[1]=boost::make_shared<Mat210x3>();
	}
}
template <int L, class real = float> Mat<L, L, real> symmetrizedMatrix(Mat<L, L, real> m) {
    size_t i, j;
    Mat<L, L, real> res = m;
    for (i = 0; i < L; ++i) {
        res[i][i] = 2 * m[i][i];
        for (j = i + 1; j < L; ++j) {
            res[i][j] = m[i][j] + m[j][i];
            res[j][i] = res[i][j];
        }
    }
    return res;
}
template< class DataTypes>
void HighOrderTriangularCorotationalFEMForceField<DataTypes>::FTCFTriangleHandler::applyCreateFunction(unsigned int triangleIndex,
        TriangleRestInformation &my_tinfo,
        const Triangle &,
        const sofa::helper::vector<unsigned int> &,
        const sofa::helper::vector<double> &)
{
	if (ff)
	{
		const std::vector< Triangle > &triangleArray=ff->_topology->getTriangles() ;
		HighOrderTriangleSetTopologyContainer *container=ff->highOrderTriangleGeo->getTopologyContainer();
		HighOrderDegreeType degree=container->getDegree();
		size_t nbControlPoints=(degree+1)*(degree+2)/2;
		size_t nbStiffnessEntries=nbControlPoints*(nbControlPoints-1)/2;
		if (my_tinfo.stiffnessVector.size()!=nbStiffnessEntries) {
			my_tinfo.stiffnessVector.resize(nbStiffnessEntries);
		}
		// set array to zero
		std::fill(my_tinfo.stiffnessVector.begin(),my_tinfo.stiffnessVector.end(),Mat2x2());

		//		const std::vector< Edge> &edgeArray=ff->_topology->getEdges() ;
		size_t i,j,k,l,m,n;

		typename DataTypes::Coord point[3];


	

		const typename DataTypes::VecCoord &restPosition=ff->mstate->read(core::ConstVecCoordId::restPosition())->getValue();
		// now computed the stiffness for the HighOrder Triangle
		sofa::helper::vector<TriangleIndexVector> tbiArray;

		tbiArray=ff->highOrderTriangleGeo->getTopologyContainer()->getTriangleIndexArray();

		size_t rank=0;
		if (ff->d_oneRotationPerIntegrationPoint.getValue()) {
			// one rotation per integration point
			Mat3x4 edgeStiffness[2];  // the off-diagonal 3x3 block matrices that makes the 12x12 linear elastic matrix
			sofa::defaulttype::Vec<3,Real> bc;
//			my_tinfo.integrationPointsStiffnessVector.resize(ff->numericalIntegrationStiffnessDataArray.size());
			my_tinfo.integrationPointsRestEdgeVector.resize(ff->numericalIntegrationStiffnessDataArray.size());
			my_tinfo.integrationPointsRestRotationArray.resize(ff->numericalIntegrationStiffnessDataArray.size());

			// loop through the integration points
			for (i=0;i<ff->numericalIntegrationStiffnessDataArray.size();++i) {

				// the barycentric coordinate
				bc=ff->numericalIntegrationStiffnessDataArray[i].integrationPoint;
				// Compute the local triangle by storing in point the 4 derivatives of the shape functions  
				ff->highOrderTriangleGeo->computeNodalValueDerivatives(triangleIndex,bc, restPosition,point);

				// initialize rotation of element
				if (ff->decompositionMethod==QR_DECOMPOSITION) {
					helper::vector<Coord> restEdgeVector(3);				
					Coord restEdgeVectorTmp[3];

					for(j=0; j<3; ++j){
						k=edgesInTriangleArray[j][0];
						l=edgesInTriangleArray[j][1];

						// store the rest edge vector
						restEdgeVector[j]=point[l]-point[k];
						restEdgeVectorTmp[j]=restEdgeVector[j];
					}
					my_tinfo.integrationPointsRestEdgeVector[i]=restEdgeVector;

					// compute the rotation matrix of the initial triangle for the QR decomposition
					computeQRRotation(my_tinfo.integrationPointsRestRotationArray[i],
						restEdgeVectorTmp);	
				} else 	if (ff->decompositionMethod==POLAR_DECOMPOSITION_MODIFIED) {
					Mat2x2 Transformation;
					Transformation[0]=point[1]-point[0];
					Transformation[1]=point[2]-point[0];
					Transformation[2]=point[3]-point[0];
					helper::Decompose<Real>::polarDecomposition( Transformation, my_tinfo.integrationPointsRestRotationArray[i] );
				}


				// compute the edge stiffness associated with that local triangle
				ff->computeTriangleStiffnessEdgeMatrix(point,edgeStiffness);
			

				my_tinfo.integrationPointsStiffnessVector.push_back(edgeStiffness[0]);
				my_tinfo.integrationPointsStiffnessVector.push_back(edgeStiffness[1]);

				my_tinfo.rotatedStiffnessVector.resize(nbControlPoints*(nbControlPoints-1)/2);
			}

		} else {
			///describe the indices of the 4 triangle vertices
			const Triangle &t= triangleArray[triangleIndex];
			//    BaseMeshTopology::EdgesInTriangle te=ff->_topology->getEdgesInTriangle(triangleIndex);


			// store the point position
			for(j=0; j<3; ++j)
				point[j]=(restPosition)[t[j]];

			if (ff->decompositionMethod==QR_DECOMPOSITION) {
				for(j=0; j<3; ++j){
					k=edgesInTriangleArray[j][0];
					l=edgesInTriangleArray[j][1];

					// store the rest edge vector
					my_tinfo.restEdgeVector[j]=point[l]-point[k];
				}
				// compute the rotation matrix of the initial triangle for the QR decomposition
				computeQRRotation(my_tinfo.restRotation,my_tinfo.restEdgeVector);
			} else 	if (ff->decompositionMethod==POLAR_DECOMPOSITION_MODIFIED) {
				Mat2x2 Transformation;
				Transformation[0]=point[1]-point[0];
				Transformation[1]=point[2]-point[0];
				helper::Decompose<Real>::polarDecomposition( Transformation, my_tinfo.restRotation );
			}
			if ((ff->integrationMethod==HighOrderTriangularCorotationalFEMForceField<DataTypes>::AFFINE_ELEMENT_INTEGRATION) || 
				((ff->d_forceAffineAssemblyForAffineElements.getValue()) && ((ff->highOrderTriangleGeo->isBezierTriangleAffine(triangleIndex,restPosition  ))))) {
				Mat3x4 edgeStiffnessVectorized[2];
				ff->computeTriangleStiffnessEdgeMatrix(point,edgeStiffnessVectorized);
				helper::system::thread::ctime_t startUpdateMat=helper::system::thread::CTime::getTime();
				if (degree==1) {
					for (rank=0;rank<nbStiffnessEntries;rank++) {
						my_tinfo.stiffnessVector[rank]+= Mat2x2((const Real *) &edgeStiffnessVectorized[0][rank][0]);
					}
				} else if (degree==2) {
					Mat15x4 res=(*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayQuadratic[0]))*edgeStiffnessVectorized[0]+
						(*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayQuadratic[1]))*edgeStiffnessVectorized[1];
					for (rank=0;rank<nbStiffnessEntries;rank++) {
						my_tinfo.stiffnessVector[rank]+= Mat2x2((const Real *) &res[rank][0]);
					}
				} else if (degree==3) {
					Mat45x4 res=(*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayCubic[0]))*edgeStiffnessVectorized[0]+
						(*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayCubic[1]))*edgeStiffnessVectorized[1];
					for (rank=0;rank<nbStiffnessEntries;rank++) {
						my_tinfo.stiffnessVector[rank]+= Mat2x2((const Real *) &res[rank][0]);
					}

				} else if (degree==4) {
					Mat105x4 res=(*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayQuartic[0]))*edgeStiffnessVectorized[0]+
						(*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayQuartic[1]))*edgeStiffnessVectorized[1];
					for (rank=0;rank<nbStiffnessEntries;rank++) {
						my_tinfo.stiffnessVector[rank]+= Mat2x2((const Real *) &res[rank][0]);
					}
				} else if (degree==5) {
					Mat210x4 res=(*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayQuintic[0]))*edgeStiffnessVectorized[0]+
						(*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayQuintic[1]))*edgeStiffnessVectorized[1];
					for (rank=0;rank<nbStiffnessEntries;rank++) {
						my_tinfo.stiffnessVector[rank]+= Mat2x2((const Real *) &res[rank][0]);
					}
				} else {
					for (rank=0;rank<nbStiffnessEntries;rank++) {
						const Vec3  & coeffVec1=ff->affineStiffnessCoefficientArray[2*rank];
						const Vec3  & coeffVec2=ff->affineStiffnessCoefficientArray[2*rank+1];

						//			Vec4 res=edgeStiffness[0]*coeffVec1+edgeStiffness[1]*coeffVec2;
						Vec4 res=edgeStiffnessVectorized[0].multTranspose(coeffVec1)+edgeStiffnessVectorized[1].multTranspose(coeffVec2);
						my_tinfo.stiffnessVector[rank]+= Mat2x2((const Real *) &res[0]);
					}
				}

				if (ff->f_printLog.getValue()) {
					helper::system::thread::ctime_t endUpdateMat=helper::system::thread::CTime::getTime();
					ff->totalUpdateMat+=endUpdateMat-startUpdateMat;
				}

			} else if (ff->integrationMethod==HighOrderTriangularCorotationalFEMForceField<DataTypes>::STANDARD_INTEGRATION){
				sofa::defaulttype::Vec<3,Real> bc;

				// loop through the integration points
				for (i=0;i<ff->numericalIntegrationStiffnessDataArray.size();++i) {

					// the barycentric coordinate
					bc=ff->numericalIntegrationStiffnessDataArray[i].integrationPoint;
					// Compute the local triangle by storing in point the 4 derivatives of the shape functions  
					ff->highOrderTriangleGeo->computeNodalValueDerivatives(triangleIndex,bc, restPosition,point);

					Mat2x2 Jacobian,inverseJacobian;
					for (j=0;j<2;++j) {
						for (k=0;k<2;++k) {
							Jacobian[j][k]=point[j][k]-point[2][k];
						}
					}
					invertMatrix(inverseJacobian,Jacobian);
					Real jac=fabs(determinant(Jacobian))*ff->numericalIntegrationStiffnessDataArray[i].integrationWeight;
					helper::system::thread::ctime_t startUpdateMat=helper::system::thread::CTime::getTime();
					helper::vector<Mat3x2> SDArray;
					for (j=0;j<nbControlPoints;j++) {
						Coord sv=inverseJacobian*ff->numericalIntegrationStiffnessDataArray[i].coefficientArray[j];
						Mat3x2 strainDisplacement;
						strainDisplacement[0][0]=sv[0];strainDisplacement[1][1]=sv[1];
						strainDisplacement[2][1]=sv[0];strainDisplacement[2][0]=sv[1];
						SDArray.push_back(strainDisplacement);	
					}

					for (rank=0,j=0;j<nbControlPoints;j++) {

						for (k=j+1;k<nbControlPoints;k++,rank++) {
							Mat2x2 edgeStiffness=(SDArray[j].transposed()*((ff->elasticityTensor)*SDArray[k]));
							edgeStiffness*= jac;
							my_tinfo.stiffnessVector[rank]+=edgeStiffness.transposed();
						}
					}
					if (ff->f_printLog.getValue()) {
						helper::system::thread::ctime_t endUpdateMat=helper::system::thread::CTime::getTime();
						ff->totalUpdateMat+=endUpdateMat-startUpdateMat;
					}
				}

			} else if (ff->integrationMethod==HighOrderTriangularCorotationalFEMForceField<DataTypes>::NUMERICAL_INTEGRATION){

				sofa::defaulttype::Vec<3,Real> bc;

				size_t p,q;
				Mat3x4 edgeStiffnessVectorized[2];
				// loop through the integration points
				for (i=0;i<ff->numericalIntegrationStiffnessDataArray.size();++i) {
					// copy weight array locally to speed-up computation
					//		std::vector<Mat4x4>  weightArray=ff->numericalIntegrationStiffnessDataArray[i].weightArray;
					// the barycentric coordinate
					bc=ff->numericalIntegrationStiffnessDataArray[i].integrationPoint;
					// Compute the local triangle by storing in point the 4 derivatives of the shape functions  
					ff->highOrderTriangleGeo->computeNodalValueDerivatives(triangleIndex,bc, restPosition,point);
					// compute the edge stiffness associated with that local triangle
					ff->computeTriangleStiffnessEdgeMatrix(point,edgeStiffnessVectorized);
					helper::system::thread::ctime_t startUpdateMat=helper::system::thread::CTime::getTime();
					// compute the stiffness matrix for all pairs of control points 
					
					if (degree==2) {
						Mat15x4 res=(*(ff->numericalIntegrationStiffnessDataArray[i].arrayPointer.weightArrayQuadratic[0]))*edgeStiffnessVectorized[0]+
							(*(ff->numericalIntegrationStiffnessDataArray[i].arrayPointer.weightArrayQuadratic[1]))*edgeStiffnessVectorized[1];
						for (rank=0;rank<nbStiffnessEntries;rank++) {
							my_tinfo.stiffnessVector[rank]+= Mat2x2((const Real *) &res[rank][0]);
						}
					} 	else if (degree==3) {
						Mat45x4 res=(*(ff->numericalIntegrationStiffnessDataArray[i].arrayPointer.weightArrayCubic[0]))*edgeStiffnessVectorized[0]+
							(*(ff->numericalIntegrationStiffnessDataArray[i].arrayPointer.weightArrayCubic[1]))*edgeStiffnessVectorized[1];
						for (rank=0;rank<nbStiffnessEntries;rank++) {
							my_tinfo.stiffnessVector[rank]+= Mat2x2((const Real *) &res[rank][0]);
						}
					} 	else if (degree==4) {
						Mat105x4 res=(*(ff->numericalIntegrationStiffnessDataArray[i].arrayPointer.weightArrayQuartic[0]))*edgeStiffnessVectorized[0]+
							(*(ff->numericalIntegrationStiffnessDataArray[i].arrayPointer.weightArrayQuartic[1]))*edgeStiffnessVectorized[1];
						for (rank=0;rank<nbStiffnessEntries;rank++) {
							my_tinfo.stiffnessVector[rank]+= Mat2x2((const Real *) &res[rank][0]);
						}
					} 	else if (degree==5) {
						Mat210x4 res=(*(ff->numericalIntegrationStiffnessDataArray[i].arrayPointer.weightArrayQuintic[0]))*edgeStiffnessVectorized[0]+
							(*(ff->numericalIntegrationStiffnessDataArray[i].arrayPointer.weightArrayQuintic[1]))*edgeStiffnessVectorized[1];
						for (rank=0;rank<nbStiffnessEntries;rank++) {
							my_tinfo.stiffnessVector[rank]+= Mat2x2((const Real *) &res[rank][0]);
						}
					} else {
						const std::vector<Vec3> & weightArray=ff->numericalIntegrationStiffnessDataArray[i].weightArray;
						for (rank=0;rank<nbStiffnessEntries;rank++) {
							const Vec3  & coeffVec1=weightArray[2*rank];
							const Vec3  & coeffVec2=weightArray[2*rank+1];

							//			Vec4 res=edgeStiffness[0]*coeffVec1+edgeStiffness[1]*coeffVec2;
							Vec4 res=edgeStiffnessVectorized[0].multTranspose(coeffVec1)+edgeStiffnessVectorized[1].multTranspose(coeffVec2);
							my_tinfo.stiffnessVector[rank]+= Mat2x2((const Real *) &res[0]);
						}
					}					

					if (ff->f_printLog.getValue()) {
						helper::system::thread::ctime_t endUpdateMat=helper::system::thread::CTime::getTime();
						ff->totalUpdateMat+=endUpdateMat-startUpdateMat;
					}

				}

			
		
//std::cerr<<" total update mat="<<((totalAssembly)/(Real)helper::system::thread::CTime::getRefTicksPerSec()) << std::endl;

			}		
            else if (ff->integrationMethod == HighOrderTriangularCorotationalFEMForceField<DataTypes>::BEZIER_NUMERICAL_INTEGRATION) {
                helper::system::thread::ctime_t startUpdateMat = helper::system::thread::CTime::getTime();
                sofa::defaulttype::Vec<3, Real> bc;
                Mat3x4 edgeStiffnessVectorized[2];
                std::vector<  sofa::defaulttype::Mat<3, 4, Real> > reducedStiffness;

                assert(ff->numericalIntegrationStiffnessDataArray.size() > 0);
                size_t numberReducedEntries = ff->numericalIntegrationStiffnessDataArray[0].weightBezierArray.size();
                reducedStiffness.resize(numberReducedEntries);

                // loop through the integration points
                for (i = 0; i < ff->numericalIntegrationStiffnessDataArray.size(); ++i) {

                    // the barycentric coordinate
                    bc = ff->numericalIntegrationStiffnessDataArray[i].integrationPoint;
                    // Compute the local tetrahedron by storing in point the 4 derivatives of the shape functions  
                    ff->highOrderTriangleGeo->computeNodalValueDerivatives(triangleIndex, bc, restPosition, point);
                    // compute the edge stiffness associated with that local tetrahedron
                    ff->computeTriangleStiffnessEdgeMatrix(point, edgeStiffnessVectorized);
                    for (j = 0; j < numberReducedEntries; ++j) {
                        reducedStiffness[j] += ff->numericalIntegrationStiffnessDataArray[i].weightBezierArray[j] * edgeStiffnessVectorized[0];
                    }
                }
                size_t r;
                for (i = 0; i < nbStiffnessEntries; ++i) {
                    sofa::defaulttype::Vec<9, int>  &mapping = ff->bezierMappingArray[i];
                    for (rank = 0, j = 0; j < 3; ++j) {
                        for (k = j + 1; k < 3; ++k, ++rank) {
                            r = j * 3 + k;
                            if (mapping[r] >= 0)
                                my_tinfo.stiffnessVector[i] += ff->bezierCoefficientArray[i][r] * Mat2x2(&reducedStiffness[(size_t)mapping[r]][rank][0]);
                            r = k * 3 + j;
                            if (mapping[r] >= 0)
                                my_tinfo.stiffnessVector[i] += ff->bezierCoefficientArray[i][r] * Mat2x2(&reducedStiffness[(size_t)mapping[r]][rank][0]).transposed();
                            r = j * 3 + j;
                            if (mapping[r] >= 0)
                                my_tinfo.stiffnessVector[i] -= 0.5*ff->bezierCoefficientArray[i][r] * symmetrizedMatrix<2, Real>(Mat2x2(&reducedStiffness[(size_t)mapping[r]][rank][0]));
                            r = k * 3 + k;
                            if (mapping[r] >= 0)
                                my_tinfo.stiffnessVector[i] -= 0.5*ff->bezierCoefficientArray[i][r] * symmetrizedMatrix<2, Real>(Mat2x2(&reducedStiffness[(size_t)mapping[r]][rank][0]));
                        }
                    }
                }
                if (ff->f_printLog.getValue()) {
                    helper::system::thread::ctime_t endUpdateMat = helper::system::thread::CTime::getTime();
                    ff->totalUpdateMat += endUpdateMat - startUpdateMat;
                }
            }
#ifdef _DEBUG
			if (ff->f_printLog.getValue()) {
				std::cerr<< " lambda="<<ff->getLambda() << " mu="<<ff->getMu() << std::endl;


				for (rank=0,l=0;l<tbiArray.size();++l)
				{
					for (m=l+1;m<tbiArray.size();++m,++rank)
					{
						std::cerr<< "Stiffness entry ["<<(unsigned int)tbiArray[l][0]<<" "<< (unsigned int)tbiArray[l][1]<<" "<< (unsigned int)tbiArray[l][2]<< "]["<<
							(unsigned int)tbiArray[m][0]<<" "<< (unsigned int)tbiArray[m][1]<<" "<< (unsigned int)tbiArray[m][2]<<"]="<<my_tinfo.stiffnessVector[rank]<<std::endl;

					}
				}
			}
#endif
		}
	}

}

template <class DataTypes> HighOrderTriangularCorotationalFEMForceField<DataTypes>::HighOrderTriangularCorotationalFEMForceField()
    : triangleInfo(initData(&triangleInfo, "triangleInfo", "Internal triangle data"))
    , _initialPoints(0)
    , updateMatrix(true)
    , d_method(initData(&d_method,std::string("linear"),"method","method for rotation computation :\"qr\" (by QR) or \"polar\" or \"polar2\" or \"none\" (Linear elastic)"))
    , d_poissonRatio(initData(&d_poissonRatio,(Real)0.3,"poissonRatio","Poisson ratio in Hooke's law"))
	, d_youngModulus(initData(&d_youngModulus,(Real)1000.,"youngModulus","Young modulus in Hooke's law"))
	, d_anisotropy(initData(&d_anisotropy,std::string("isotropic"),"elasticitySymmetry","the type of anisotropy for the elasticity tensor :\"isotropic\"  or \"transverseIsotropic\" or \"orthotropic\" or \"cubic\" "))
	, d_anisotropyParameter(initData(&d_anisotropyParameter,"anisotropyParameters","the elastic parameters for anisotropic materials "))
	, d_anisotropyDirection(initData(&d_anisotropyDirection,"anisotropyDirections","the directions of anisotropy"))
	, numericalIntegrationOrder( initData(&numericalIntegrationOrder,(size_t)2,"integrationOrder","The order of integration for numerical integration"))
	, d_integrationMethod( initData(&d_integrationMethod,std::string("analytical"),"integrationMethod","\"analytical\" if closed form expression for affine element, \"numerical\" if numerical integration is chosen,  \"standard\" if standard integration is chosen"))
	, numericalIntegrationMethod( initData(&numericalIntegrationMethod, std::string("Triangle Gauss"),"numericalIntegrationMethod","The type of numerical integration method chosen"))
	 , d_oneRotationPerIntegrationPoint(initData(&d_oneRotationPerIntegrationPoint,false,"oneRotationPerIntegrationPoint","if true then computes one rotation per integration point"))
	 , d_assemblyTime(initData(&d_assemblyTime,(Real)0,"assemblyTime","the time spent in assembling the stiffness matrix. Only updated if printLog is set to true"))
	 , d_forceAffineAssemblyForAffineElements(initData(&d_forceAffineAssemblyForAffineElements,true,"forceAffineAssemblyForAffineElements","if true affine triangles are always assembled with the closed form formula, Otherwise use the method defined in integrationMethod"))
	 , lambda(0)
    , mu(0)
    , triangleHandler(NULL)
{
    triangleHandler = new FTCFTriangleHandler(this,&triangleInfo);
}

template <class DataTypes> HighOrderTriangularCorotationalFEMForceField<DataTypes>::~HighOrderTriangularCorotationalFEMForceField()
{
    if (triangleHandler) delete triangleHandler;
}

 template <class DataTypes> void HighOrderTriangularCorotationalFEMForceField<DataTypes>::assembleAnisotropicTensors()
 {

 }
template <class DataTypes> void HighOrderTriangularCorotationalFEMForceField<DataTypes>::init()
{
    //	serr << "initializing HighOrderTriangularCorotationalFEMForceField" << sendl;
    this->Inherited::init();

    _topology = this->getContext()->getMeshTopology();
	this->getContext()->get(highOrderTriangleGeo);

    if ((_topology->getNbTriangles()==0) || (!highOrderTriangleGeo))
    {
        serr << "ERROR(HighOrderTriangularCorotationalFEMForceField): object must have a Triangular Set Topology and a HighOrderTriangleSetGeometryAlgorithms component "<<sendl;
        return;
    }
    updateLameCoefficients();

	if (d_anisotropy.getValue() == "isotropic")
        elasticitySymmetry= ISOTROPIC;
    else if (d_anisotropy.getValue() == "cubic") 
        elasticitySymmetry= CUBIC;
    else if (d_anisotropy.getValue() == "transverseIsotropic") 
        elasticitySymmetry= TRANSVERSE_ISOTROPIC;
	else
    {
        serr << "cannot recognize symmetry "<< d_anisotropy.getValue() << ". Must be either \"isotropic\" or \"cubic\"  or \"transverseIsotropic\" or \"orthotropic\"" << sendl;
    }

    if (d_method.getValue() == "polar")
        decompositionMethod= POLAR_DECOMPOSITION;
    else if ((d_method.getValue() == "qr") || (d_method.getValue() == "large"))
        decompositionMethod= QR_DECOMPOSITION;
    else if (d_method.getValue() == "polar2")
        decompositionMethod= POLAR_DECOMPOSITION_MODIFIED;
	 else if ((d_method.getValue() == "none") || (d_method.getValue() == "linear"))
        decompositionMethod= LINEAR_ELASTIC;
    else
    {
        serr << "cannot recognize method "<< d_method.getValue() << ". Must be either qr (or large), polar, polar2 or none (or linear)" << sendl;
    }
    if (d_integrationMethod.getValue() == "analytical")
        integrationMethod = AFFINE_ELEMENT_INTEGRATION;
    else if (d_integrationMethod.getValue() == "numerical")
        integrationMethod = NUMERICAL_INTEGRATION;
    else if (d_integrationMethod.getValue() == "numerical2")
        integrationMethod = NUMERICAL_INTEGRATION_2;
    else if (d_integrationMethod.getValue() == "bezierNumerical")
        integrationMethod = BEZIER_NUMERICAL_INTEGRATION;
    else if (d_integrationMethod.getValue() == "standard")
        integrationMethod = STANDARD_INTEGRATION;
    else
    {
        serr << "cannot recognize method " << d_integrationMethod.getValue() << ". Must be either \"analytical\" or \"numerical\" or \"bezierNumerical\" or \"standard\"" << sendl;
    }
    std::set<typename topology::NumericalIntegrationDescriptor<Real, 3>::QuadratureMethod> qmSet = highOrderTriangleGeo->getTriangleNumericalIntegrationDescriptor().getQuadratureMethods();
    if (qmSet.count(numericalIntegrationMethod.getValue()) == 0) {
        serr << "cannot recognize numerical integration method  " << numericalIntegrationMethod.getValue() << sendl;
    }
    helper::vector<TriangleRestInformation>& triangleInf = *(triangleInfo.beginEdit());
    triangleInf.resize(_topology->getNbTriangles());

    if (_initialPoints.size() == 0)
    {
        // get restPosition
        const VecCoord& p = this->mstate->read(core::ConstVecCoordId::restPosition())->getValue();
        _initialPoints=p;
	}

	helper::system::thread::ctime_t startAssembly=helper::system::thread::CTime::getTime();
	totalUpdateMat=0;

    size_t i;


	// precompute the coefficients for handling affine elements
	topology::TriangleIndexVector tbi1,tbi2;
	affineStiffnessCoefficientArray.clear();
	
	topology::HighOrderDegreeType degree=highOrderTriangleGeo->getTopologyContainer()->getDegree();
	size_t nbControlPoints=(degree+1)*(degree+2)/2;
	sofa::helper::vector<topology::TriangleIndexVector> tbiArray;
	tbiArray=highOrderTriangleGeo->getTopologyContainer()->getTriangleIndexArray();

	if (degree==1) 
		integrationMethod= AFFINE_ELEMENT_INTEGRATION;

	computeElasticityTensor();
	if (degree>1) {
		if ((integrationMethod== AFFINE_ELEMENT_INTEGRATION) || (d_forceAffineAssemblyForAffineElements.getValue()))
		{
			if (degree<6) {
				affineStiffnessCoefficientPreStoredArray.allocate(degree);
			}
			std::vector<Real> coeffArray(3);
			Mat3x3 coeffMatrix;
			size_t j,k,l,m,n,rank;
			for (rank=0,j=0;j<nbControlPoints;j++) {
				tbi1=tbiArray[j];
				for (k=j+1;k<nbControlPoints;k++,rank++) {
					tbi2=tbiArray[k];
					coeffMatrix=highOrderTriangleGeo->getAffineStiffnessCoefficientMatrix(tbi1,tbi2);
					// substract the diagonal terms such that only edge stiffness are used
					for(l=0; l<3; ++l){
						for(m=0; m<3; ++m){
							if (m!=l) {
								coeffMatrix[l][m]-=0.5*(coeffMatrix[l][l]+coeffMatrix[m][m]);
							}
						}
					}
					if (degree==2) {
						for(l=0; l<3; ++l){
							m=edgesInTriangleArray[l][0];
							n=edgesInTriangleArray[l][1];

							(*(affineStiffnessCoefficientPreStoredArray.weightArrayQuadratic[0]))[rank][l]=coeffMatrix[m][n];
							(*(affineStiffnessCoefficientPreStoredArray.weightArrayQuadratic[1]))[rank][l]=coeffMatrix[n][m];
						}
					} else if (degree==3) {
						for(l=0; l<3; ++l){
							m=edgesInTriangleArray[l][0];
							n=edgesInTriangleArray[l][1];

							(*(affineStiffnessCoefficientPreStoredArray.weightArrayCubic[0]))[rank][l]=coeffMatrix[m][n];
							(*(affineStiffnessCoefficientPreStoredArray.weightArrayCubic[1]))[rank][l]=coeffMatrix[n][m];
						}
					} else if (degree==4) {
						for(l=0; l<3; ++l){
							m=edgesInTriangleArray[l][0];
							n=edgesInTriangleArray[l][1];

							(*(affineStiffnessCoefficientPreStoredArray.weightArrayQuartic[0]))[rank][l]=coeffMatrix[m][n];
							(*(affineStiffnessCoefficientPreStoredArray.weightArrayQuartic[1]))[rank][l]=coeffMatrix[n][m];
						}
					} else if (degree==5) {
						for(l=0; l<3; ++l){
							m=edgesInTriangleArray[l][0];
							n=edgesInTriangleArray[l][1];

							(*(affineStiffnessCoefficientPreStoredArray.weightArrayQuintic[0]))[rank][l]=coeffMatrix[m][n];
							(*(affineStiffnessCoefficientPreStoredArray.weightArrayQuintic[1]))[rank][l]=coeffMatrix[n][m];
						}

					} else {
						Vec3 coeffVec1,coeffVec2;
						for(l=0; l<3; ++l){
							m=edgesInTriangleArray[l][0];
							n=edgesInTriangleArray[l][1];
							coeffVec1[l]=coeffMatrix[m][n];
							coeffVec2[l]=coeffMatrix[n][m];
						}

						affineStiffnessCoefficientArray.push_back(coeffVec1);
						affineStiffnessCoefficientArray.push_back(coeffVec2);
					}
				}
			}
		}
		if (integrationMethod== NUMERICAL_INTEGRATION)  
		{
			numericalIntegrationStiffnessDataArray.clear();
			/// get value of integration points0
			topology::NumericalIntegrationDescriptor<Real,3> &nid=highOrderTriangleGeo->getTriangleNumericalIntegrationDescriptor();
			typename topology::NumericalIntegrationDescriptor<Real,3>::QuadraturePointArray qpa=nid.getQuadratureMethod((typename topology::NumericalIntegrationDescriptor<Real,3>::QuadratureMethod)numericalIntegrationMethod.getValue(),
				numericalIntegrationOrder.getValue());
			size_t i,j,k,l,m,n;
			sofa::defaulttype::Vec<3,Real> bc;
			Real weight;
			Mat3x3 coeffMatrix;

			// loop through the integration points
			for (i=0;i<qpa.size();++i) {
				NumericalIntegrationStiffnessData nimd;
				typename topology::NumericalIntegrationDescriptor<Real,3>::QuadraturePoint qp=qpa[i];
				// the barycentric coordinate
				nimd.integrationPoint=qp.first;
				// the weight of the integration point
				weight=qp.second;
				nimd.integrationWeight=qp.second;
				if ((integrationMethod== NUMERICAL_INTEGRATION) && (degree<6)) {
					nimd.arrayPointer.allocate(degree);
				}
				std::vector<Vec3> shapeFunctionDerivativeArray;
				for(j=0;j<tbiArray.size();++j) {
					Vec3 deriv=highOrderTriangleGeo->computeShapeFunctionDerivatives(tbiArray[j],qp.first);
					shapeFunctionDerivativeArray.push_back(deriv);
					Deriv der(deriv[0]-deriv[2],deriv[1]-deriv[2]);
					nimd.coefficientArray.push_back(der);
				}
				size_t rank;
				for(rank=0,j=0;j<tbiArray.size();++j) {
					for(k=j+1;k<tbiArray.size();++k,++rank) {
						coeffMatrix=dyad(shapeFunctionDerivativeArray[j],shapeFunctionDerivativeArray[k])*2*weight;
						for(l=0; l<3; ++l){
							for(m=0; m<3; ++m){
								if (l!=m) {
									coeffMatrix[l][m]-=0.5*(coeffMatrix[l][l]+coeffMatrix[m][m]);
								}
							}
						}
						if (integrationMethod== NUMERICAL_INTEGRATION) { 
							if ((degree>5) || (d_oneRotationPerIntegrationPoint.getValue())) {
								Vec3 coeffVec1,coeffVec2;
								for(l=0; l<3; ++l){
									m=edgesInTriangleArray[l][0];
									n=edgesInTriangleArray[l][1];
									coeffVec1[l]=coeffMatrix[m][n];
									coeffVec2[l]=coeffMatrix[n][m];
								}
								nimd.weightArray.push_back(coeffVec1);
								nimd.weightArray.push_back(coeffVec2);
							} else {
								if (degree==2) {

									for(l=0; l<3; ++l){
										m=edgesInTriangleArray[l][0];
										n=edgesInTriangleArray[l][1];

										(*(nimd.arrayPointer.weightArrayQuadratic[0]))[rank][l]=coeffMatrix[m][n];
										(*(nimd.arrayPointer.weightArrayQuadratic[1]))[rank][l]=coeffMatrix[n][m];
									}
								} else 	if (degree==3) {

									for(l=0; l<3; ++l) {
										m=edgesInTriangleArray[l][0];
										n=edgesInTriangleArray[l][1];

										(*(nimd.arrayPointer.weightArrayCubic[0]))[rank][l]=coeffMatrix[m][n];
										(*(nimd.arrayPointer.weightArrayCubic[1]))[rank][l]=coeffMatrix[n][m];
									}
								} else 	if (degree==4) {

									for(l=0; l<3; ++l) {
										m=edgesInTriangleArray[l][0];
										n=edgesInTriangleArray[l][1];

										(*(nimd.arrayPointer.weightArrayQuartic[0]))[rank][l]=coeffMatrix[m][n];
										(*(nimd.arrayPointer.weightArrayQuartic[1]))[rank][l]=coeffMatrix[n][m];
									}
								} else 	if (degree==5) {
									for(l=0; l<3; ++l) {
										m=edgesInTriangleArray[l][0];
										n=edgesInTriangleArray[l][1];

										(*(nimd.arrayPointer.weightArrayQuintic[0]))[rank][l]=coeffMatrix[m][n];
										(*(nimd.arrayPointer.weightArrayQuintic[1]))[rank][l]=coeffMatrix[n][m];
									}
								}
							}
						}
					}
				}

				numericalIntegrationStiffnessDataArray.push_back(nimd);
			}
		}
        if (integrationMethod == BEZIER_NUMERICAL_INTEGRATION)
        {
            /// first fill the first vector independent from the integration points
            BezierTriangleSetGeometryAlgorithms<DataTypes> *bezierTriangleGeo = dynamic_cast<BezierTriangleSetGeometryAlgorithms<DataTypes> *> (highOrderTriangleGeo);
            if (bezierTriangleGeo == NULL) {
                serr << "Could not find any BezierTriangleSetGeometryAlgorithms while using BEZIER_NUMERICAL_INTEGRATION" << sendl;
                return;
            }
            bezierCoefficientArray.clear();
            /// store the coefficient that are independent from the integration point
            topology::TetrahedronIndexVector tbi1Copy, tbi2Copy;
            size_t rank, r;
            size_t i, j, k, l, m;
            for (rank = 0, j = 0; j < tbiArray.size(); ++j) {
                for (k = j + 1; k < tbiArray.size(); ++k, ++rank) {
                    Vec9 weight;
                    tbi1 = tbiArray[j];
                    tbi2 = tbiArray[k];
                    for (r = 0, l = 0; l < 3; ++l) {
                        for (m = 0; m < 3; ++m, ++r) {
                            if ((tbi1[l] * tbi2[m]) != 0) {
                                tbi1Copy = tbi1;
                                tbi1Copy[l] -= 1;
                                weight[r] = 1.0f / (factorialTVI(tbi1Copy));
                                tbi2Copy = tbi2;
                                tbi2Copy[m] -= 1;
                                weight[r] *= 1.0f / (factorialTVI(tbi2Copy));
                            }
                            else
                                weight[r] = 0;

                        }
                    }
                    bezierCoefficientArray.push_back(weight);
                }
            }
            /// store the mapping coefficient that are independent from the integration point
            // now fills the vector for each integration point
            sofa::helper::vector<topology::TriangleIndexVector> tbiArray2;
            tbiArray2 = bezierTriangleGeo->getTopologyContainer()->getTriangleIndexArrayOfGivenDegree(2 * degree - 2);

            // first create a map to speed up the assignment of index from  TetrahedronIndexVector
            std::map<topology::TriangleIndexVector, size_t> tivMap;
            std::map<topology::TriangleIndexVector, size_t>::iterator itmap;
            for (j = 0; j < tbiArray2.size(); ++j) {
                tivMap.insert(std::make_pair(tbiArray2[j], j));
            }
            /// now fills the bezierMappingArray
            for (rank = 0, j = 0; j < tbiArray.size(); ++j) {
                for (k = j + 1; k < tbiArray.size(); ++k, ++rank) {
                    tbi1 = tbiArray[j] + tbiArray[k];
                    Vec9Int mapping;
                    for (r = 0, l = 0; l<3; ++l) {
                        for (m = 0; m<3; ++m, ++r) {

                            if (((tbi1[l] * tbi1[m]) == 0) || (((l == m) && (tbi1[l] < 2)))) {
                                //        if ((tbi1[l] * tbi1[m]) == 0) {
                                mapping[r] = -1;
                            }
                            else {

                                tbi2 = tbi1;
                                tbi2[l] -= 1;
                                tbi2[m] -= 1;
                                itmap = tivMap.find(tbi2);
                                assert(itmap != tivMap.end());
                                mapping[r] = (*itmap).second;

                            }
                        }
                    }
                    bezierMappingArray.push_back(mapping);
                }
            }

            numericalIntegrationStiffnessDataArray.clear();
            /// get value of integration points0
            topology::NumericalIntegrationDescriptor<Real, 3> &nid = highOrderTriangleGeo->getTriangleNumericalIntegrationDescriptor();
            typename topology::NumericalIntegrationDescriptor<Real, 3>::QuadraturePointArray qpa = nid.getQuadratureMethod((typename topology::NumericalIntegrationDescriptor<Real, 3>::QuadratureMethod)numericalIntegrationMethod.getValue(),
                numericalIntegrationOrder.getValue());

            sofa::defaulttype::Vec<4, Real> bc;
            Real weight, fac;
           
            fac = (Real)lfactorial(degree - 1)*(Real)lfactorial(degree - 1) / (Real)lfactorial(2 * degree - 2);

            // loop through the integration points
            for (i = 0; i<qpa.size(); ++i) {
                NumericalIntegrationStiffnessData nimd;
                typename topology::NumericalIntegrationDescriptor<Real, 3>::QuadraturePoint qp = qpa[i];
                // the barycentric coordinate
                nimd.integrationPoint = qp.first;
                // the weight of the integration point
                weight = qp.second;
                nimd.integrationWeight = qp.second;

                nimd.weightBezierArray.resize(tbiArray2.size());
                for (j = 0; j < tbiArray2.size(); ++j) {
                    nimd.weightBezierArray[j] = 2 * degree*degree*fac*bezierTriangleGeo->computeShapeFunctionOfGivenDegree(tbiArray2[j], qp.first, 2 * degree - 2);
                    nimd.weightBezierArray[j] *= weight*factorialTVI(tbiArray2[j]);
                }

                numericalIntegrationStiffnessDataArray.push_back(nimd);
            }
        }
		if (integrationMethod== STANDARD_INTEGRATION) 
		{
			numericalIntegrationStiffnessDataArray.clear();
			/// get value of integration points0
			topology::NumericalIntegrationDescriptor<Real,3> &nid=highOrderTriangleGeo->getTriangleNumericalIntegrationDescriptor();
			typename topology::NumericalIntegrationDescriptor<Real,3>::QuadraturePointArray qpa=nid.getQuadratureMethod((typename topology::NumericalIntegrationDescriptor<Real,3>::QuadratureMethod)numericalIntegrationMethod.getValue(),
				numericalIntegrationOrder.getValue());
			size_t i,j,k;
			sofa::defaulttype::Vec<3,Real> bc;
			Real weight;
			Mat3x3 coeffMatrix;

			// loop through the integration points
			for (i=0;i<qpa.size();++i) {
				NumericalIntegrationStiffnessData nimd;
				typename topology::NumericalIntegrationDescriptor<Real,3>::QuadraturePoint qp=qpa[i];
				// the barycentric coordinate
				nimd.integrationPoint=qp.first;
				// the weight of the integration point	
				nimd.integrationWeight=qp.second;

				std::vector<Vec3> shapeFunctionDerivativeArray;
				for(j=0;j<tbiArray.size();++j) {
					Vec3 deriv=highOrderTriangleGeo->computeShapeFunctionDerivatives(tbiArray[j],qp.first);
					Deriv der(deriv[0]-deriv[2],deriv[1]-deriv[2]);
					nimd.coefficientArray.push_back(der);
				}
				numericalIntegrationStiffnessDataArray.push_back(nimd);
			}



		}
	}
	    /// initialize the data structure associated with each triangle
    for (i=0; i<_topology->getNbTriangles(); ++i)
    {
        triangleHandler->applyCreateFunction(i,triangleInf[i],_topology->getTriangle(i),
                (const helper::vector< unsigned int > )0,
                (const helper::vector< double >)0);
    }
	if (this->f_printLog.getValue()) {
		helper::system::thread::ctime_t endAssembly=helper::system::thread::CTime::getTime();
		std::cerr<< "Assembly time="<< ((endAssembly-startAssembly)/(Real)helper::system::thread::CTime::getRefTicksPerSec()) << std::endl;
		std::cerr<<" total update mat="<<((totalUpdateMat)/(Real)helper::system::thread::CTime::getRefTicksPerSec()) << std::endl;
		d_assemblyTime.setValue(((endAssembly-startAssembly)/(Real)helper::system::thread::CTime::getRefTicksPerSec()));
	}
    /// set the call back function upon creation of a triangle
    triangleInfo.createTopologicalEngine(_topology,triangleHandler);
    triangleInfo.registerTopologicalData();
    triangleInfo.endEdit();

	updateTopologyInfo=true;

}


template <class DataTypes>
void HighOrderTriangularCorotationalFEMForceField<DataTypes>::updateTopologyInformation()
{
    int i;
    unsigned int j;

    int nbTriangles=_topology->getNbTriangles();

    TriangleRestInformation *tetinfo;

    helper::vector<typename HighOrderTriangularCorotationalFEMForceField<DataTypes>::TriangleRestInformation>& triangleInf = *(triangleInfo.beginEdit());
   


    for(i=0; i<nbTriangles; i++ )
    {
        tetinfo=&triangleInf[i];
        /// describe the jth edge index of triangle no i
        const EdgesInTriangle &tea= _topology->getEdgesInTriangle(i);
        /// describe the jth vertex index of triangle no i
        const Triangle &ta= _topology->getTriangle(i);

        for (j=0; j<3; ++j)
        {
            tetinfo->v[j]=ta[j];
        }


    }
    updateTopologyInfo=false;
    triangleInfo.endEdit();
}
template<class DataTypes>
void HighOrderTriangularCorotationalFEMForceField<DataTypes>::computeElasticityTensor() 											
{
	const helper::vector<Real> & anisotropyParameter=d_anisotropyParameter.getValue();
	const helper::vector<Coord> & anisotropyDirection=d_anisotropyDirection.getValue();
	if (elasticitySymmetry==ISOTROPIC) {
			// elasticity tensor in isotropic case
		elasticityTensor.clear();
		Real lambda=getLambda();
		Real mu=getMu();
		elasticityTensor(0,0)=mu+lambda;elasticityTensor(1,1)=mu+lambda;
		elasticityTensor(0,1)=lambda;elasticityTensor(1,0)=lambda;
		elasticityTensor(2,2)=mu/2;
	
	} /* else if (elasticitySymmetry==CUBIC) {
		assert(anisotropyParameter.size()>=1);
		assert(anisotropyDirection.size()>=1);
		// get 3 orthogonal direction starting from the direction of anisotropy
		Coord n=anisotropyDirection[0];
		
		n/=n.norm();
		Coord v1,v2;
		if ((n[0]!=0) || (n[1]!=0)) {
			v1=Coord(-n[1],n[0],n[2]);
		} else {
			v1=Coord(1,0,0);
		}
		v1=cross(n,v1);
		v1/=v1.norm();
		v2=cross(v1,n);

		// build the orthogonal matrices
		Mat2x2 Nd;
		Nd.identity();
		Nd/= sqrt(3.0f);
		Mat2x2 Ne=(2*dyad(n,n)-dyad(v1,v1)-dyad(v2,v2))/sqrt(6.0f);
		Mat2x2 Np=(dyad(v1,v1)-dyad(v2,v2))/sqrt(2.0f);
		Mat2x2 Ns1=(dyad(v2,n)+dyad(n,v2))/sqrt(2.0f);
		Mat2x2 Ns2=(dyad(v1,n)+dyad(n,v1))/sqrt(2.0f);
		Mat2x2 Ns3=(dyad(v1,v2)+dyad(v2,v1))/sqrt(2.0f);
		// get the constants from the young modulus, Poisson ratio and anisotropy ratio.
		Real youngModulus=d_youngModulus.getValue();
		Real poissonRatio=d_poissonRatio.getValue();
		Real anisotropyRatio=anisotropyParameter[0];
		// use equations from http://solidmechanics.org/text/Chapter3_2/Chapter3_2.htm
		Real c11,c12;
		c11=youngModulus*(1-poissonRatio)/(1-poissonRatio-poissonRatio*poissonRatio);
		c12=youngModulus*poissonRatio/(1-poissonRatio-poissonRatio*poissonRatio);
		Real c44=anisotropyRatio*(c11-c12)/2;
		// push all symmetric matrices and the eigenvalues
		anisotropyMatrixArray.push_back(Nd);
		anisotropyScalarArray.push_back(c11+2*c12);
		anisotropyMatrixArray.push_back(Ne);
		anisotropyScalarArray.push_back(c11-c12);
		anisotropyMatrixArray.push_back(Np);
		anisotropyScalarArray.push_back(c11-c12);
		anisotropyMatrixArray.push_back(Ns1);
		anisotropyScalarArray.push_back(c44);
		anisotropyMatrixArray.push_back(Ns2);
		anisotropyScalarArray.push_back(c44);
		anisotropyMatrixArray.push_back(Ns3);
		anisotropyScalarArray.push_back(c44);
	} else if (elasticitySymmetry==TRANSVERSE_ISOTROPIC) {
		assert(anisotropyParameter.size()>=1);
		assert(anisotropyDirection.size()>=1);
		// get 3 orthogonal direction starting from the direction of anisotropy
		Coord n=anisotropyDirection[0];

			n/=n.norm();
		Coord v1,v2;
		if ((n[0]!=0) || (n[1]!=0)) {
			v1=Coord(-n[1],n[0],n[2]);
		} else {
			v1=Coord(1,0,0);
		}
		v1=cross(n,v1);
		v1/=v1.norm();
		v2=cross(v1,n);

		// get the constants from the young modulus, Poisson ratio and anisotropy ratio.
		Real youngModulusTransverse=d_youngModulus.getValue();
		Real poissonRatioTransverse=d_poissonRatio.getValue();
		Real youngModulusLongitudinal=anisotropyParameter[0];
		Real poissonRatioTransverseLongitudinal=anisotropyParameter[1];
		Real shearModulusTransverse=anisotropyParameter[2];
		// use equations from http://solidmechanics.org/text/Chapter3_2/Chapter3_2.htm
		Real poissonRatioLongitudinalTransverse=poissonRatioTransverseLongitudinal*youngModulusLongitudinal/youngModulusTransverse;
		Real gamma=1/(1-poissonRatioTransverse*poissonRatioTransverse-2*poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal-2*poissonRatioTransverse*poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal);
		
		Real c11=youngModulusTransverse*gamma*(1-poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal);
		Real c33=youngModulusLongitudinal*gamma*(1-poissonRatioTransverse*poissonRatioTransverse);
		Real c12=youngModulusTransverse*gamma*(poissonRatioTransverse+poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal);
		Real c13=youngModulusTransverse*gamma*(poissonRatioLongitudinalTransverse+poissonRatioTransverse*poissonRatioTransverseLongitudinal);
		Real c66=youngModulusTransverse/(2*(1+poissonRatioTransverse));
		Real c44=shearModulusTransverse;

		Real talpha=sqrt(2.0f)*(c11+c12-c33)/(4*c13);
		Real alpha=atan(talpha);
		Real salpha=sin(alpha); 
		Real calpha=cos(alpha);
		Real secalpha=1/calpha;

		Real eigen1=c33+M_SQRT2*c13*(talpha+secalpha);
		Real eigen2=c11-c12; 
		Real eigen3=c33+M_SQRT2*c13*(talpha-secalpha);
		Real eigen4=shearModulusTransverse; 

		// build the orthogonal matrices
		
		Mat2x2 tmp;
		Real val1=(0.5*(1+salpha)+sqrt(2.0)*calpha/4.0f); 
		Real val2=(0.5*(1-salpha)+sqrt(2.0)*calpha/2.0f); 
		Mat2x2 Nh1=val1*(dyad(v1,v1)+dyad(v2,v2))+val2*dyad(n,n);
		Nh1/=sqrt(2*val1*val1+val2*val2);
		val1=(0.5*(1-salpha)-sqrt(2.0)*calpha/4.0f); 
		val2=(0.5*(1+salpha)-sqrt(2.0)*calpha/2.0f); 
		Mat2x2 Nh2=val1*(dyad(v1,v1)+dyad(v2,v2))+val2*dyad(n,n);
		Nh2/=sqrt(2*val1*val1+val2*val2);

		Mat2x2 Np=(dyad(v1,v1)-dyad(v2,v2))/sqrt(2.0f);
		Mat2x2 Ns1=(dyad(v2,n)+dyad(n,v2))/sqrt(2.0f);
		Mat2x2 Ns2=(dyad(v1,n)+dyad(n,v1))/sqrt(2.0f);
		Mat2x2 Ns3=(dyad(v1,v2)+dyad(v2,v1))/sqrt(2.0f);
	
		
		
		// push all symmetric matrices and the eigenvalues
		anisotropyMatrixArray.push_back(Nh1);
		anisotropyScalarArray.push_back(eigen1);
		anisotropyMatrixArray.push_back(Np);
		anisotropyScalarArray.push_back(eigen2);
		anisotropyMatrixArray.push_back(Ns3);
		anisotropyScalarArray.push_back(eigen2);
		anisotropyMatrixArray.push_back(Nh2);
		anisotropyScalarArray.push_back(eigen3);
		anisotropyMatrixArray.push_back(Ns1);
		anisotropyScalarArray.push_back(eigen4);
		anisotropyMatrixArray.push_back(Ns2);
		anisotropyScalarArray.push_back(eigen4);
	} */
} 

template<class DataTypes>
void HighOrderTriangularCorotationalFEMForceField<DataTypes>::computeTriangleStiffnessEdgeMatrix(const Coord point[3],
																								  Mat3x4 edgeStiffnessVectorized[2])
{
	Coord shapeVector[3];
	Mat2x2 edgeStiffness[3];
	/// compute 6 times the rest volume
	Real volume=areaProduct(point[1]-point[0],point[2]-point[0]);
	/// store the rest volume
//	my_tinfo.restVolume=volume/6;

	size_t j,k,l,m,n;
	// store shape vectors at the rest configuration
	for(j=0; j<3; ++j)
	{
		shapeVector[j]=ortho(point[(j+2)%3]-point[(j+1)%3])/volume;		
	}
	if (elasticitySymmetry==ISOTROPIC) {
		Real mu=getMu()*fabs(volume)/4;
		Real lambda=getLambda()*fabs(volume)/2;
		Real val;

		/// compute the edge stiffness of the linear elastic material
		for(j=0; j<3; ++j)
		{
			k=edgesInTriangleArray[j][0];
			l=edgesInTriangleArray[j][1];
			// the linear stiffness matrix using shape vectors and Lame coefficients
			val=mu*dot(shapeVector[l],shapeVector[k]);
			for(m=0; m<2; ++m)
			{
				for(n=0; n<2; ++n)
				{
					edgeStiffness[j][m][n]=lambda*shapeVector[k][n]*shapeVector[l][m]+
						mu*shapeVector[l][n]*shapeVector[k][m];

					if (m==n)
					{
						edgeStiffness[j][m][m]+=(Real)val;
					}
				}
			}
		}
	} else {
		size_t i;
		for(j=0; j<3; ++j)
		{
			k=edgesInTriangleArray[j][0];
			l=edgesInTriangleArray[j][1];
			// the linear stiffness matrix using shape vectors and Lame coefficients
			Mat2x2 tmp=dyad(shapeVector[l],shapeVector[k]);
			for(i=0;i<anisotropyScalarArray.size();++i) {
				edgeStiffness[j]+=anisotropyScalarArray[i]*anisotropyMatrixArray[j]*tmp*anisotropyMatrixArray[j];
			}
		}
	}
	size_t p;
	for(j=0; j<3; ++j)
	{
		k=edgesInTriangleArray[j][0];
		l=edgesInTriangleArray[j][1];
		for(p=0,m=0; m<2; ++m)
		{
			for(n=0; n<2; ++n,++p)
			{
				edgeStiffnessVectorized[0][j][p]=edgeStiffness[j][m][n];
				edgeStiffnessVectorized[1][j][p]=edgeStiffness[j][n][m];
			}
		}
	}
	
}

template<class DataTypes>
void HighOrderTriangularCorotationalFEMForceField<DataTypes>::computeQRRotation( Mat2x2 &r, const Coord *dp)
{
    // first vector on first edge
    // second vector in the plane of the two first edges
    // third vector orthogonal to first and second

    Coord edgex = dp[0];
    edgex.normalize();
	/*
    Coord edgey = dp[1];

    Coord edgez = cross( edgex, edgey );
    edgez.normalize();

    edgey = cross( edgez, edgex );
    edgey.normalize();

    r[0][0] = edgex[0];
    r[0][1] = edgex[1];
    r[0][2] = edgex[2];
    r[1][0] = edgey[0];
    r[1][1] = edgey[1];
    r[1][2] = edgey[2];
    r[2][0] = edgez[0];
    r[2][1] = edgez[1];
    r[2][2] = edgez[2];*/
}

template <class DataTypes>
const helper::vector<typename HighOrderTriangularCorotationalFEMForceField<DataTypes>::Mat2x2> & 
	HighOrderTriangularCorotationalFEMForceField<DataTypes>::getStiffnessArray(
	const size_t i,
	typename HighOrderTriangularCorotationalFEMForceField<DataTypes>::TriangleRestInformation *restTriangle)
{
	return(restTriangle->stiffnessVector);
}
template <class DataTypes>
const  helper::vector<typename HighOrderTriangularCorotationalFEMForceField<DataTypes>::Mat2x2> & 
	HighOrderTriangularCorotationalFEMForceField<DataTypes>::getRotatedStiffnessArray(
	const size_t i,
	const typename HighOrderTriangularCorotationalFEMForceField<DataTypes>::TriangleRestInformation *restTriangle)
{

	if (decompositionMethod==LINEAR_ELASTIC)
	{
		return(restTriangle->stiffnessVector);
	} else
		return(restTriangle->rotatedStiffnessVector);
}
template <class DataTypes>
void HighOrderTriangularCorotationalFEMForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* mparams /* PARAMS FIRST */, 
																	 DataVecDeriv &  dataF, const DataVecCoord &  dataX , const DataVecDeriv & /*dataV*/ )
{

    VecDeriv& f        = *(dataF.beginEdit());
    const VecCoord& x  =   dataX.getValue()  ;
    const VecCoord& x0= this->mstate->read(core::ConstVecCoordId::restPosition())->getValue();
	sofa::helper::vector<Coord> dp,force;
    size_t i,j,k,l,v0,v1,rank,p;
    size_t nbTriangles=_topology->getNbTriangles();
    
	HighOrderDegreeType degree=highOrderTriangleGeo->getTopologyContainer()->getDegree();
	size_t nbControlPoints=(degree+1)*(degree+2)/2;
	HighOrderTriangleSetTopologyContainer::VecPointID indexArray;


    if (updateTopologyInfo)
    {
        updateTopologyInformation();
    }
    helper::vector<TriangleRestInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleRestInformation *tetinfo;
	
	dp.resize(nbControlPoints);
	force.resize(nbControlPoints);
	
	Coord dpos,sv;
		
	for(i=0; i<nbTriangles; i++ )
	{
		tetinfo=&triangleInf[i];
		const HighOrderTriangleSetTopologyContainer::VecPointID &indexArray=
			highOrderTriangleGeo->getTopologyContainer()->getGlobalIndexArrayOfControlPoints(i);
		
		nbControlPoints=indexArray.size();
	
		if (d_oneRotationPerIntegrationPoint.getValue()) {
			
			Mat2x2 S,R;
			size_t j,k,l,m,n;
//			helper::vector<Mat2x2> stiffnessArray(nbControlPoints*(nbControlPoints-1)/2);
			helper::vector<Mat2x2> &stiffnessArray=tetinfo->rotatedStiffnessVector;
			std::fill(stiffnessArray.begin(),stiffnessArray.end(),Mat2x2());
			assert(stiffnessArray.size()==nbControlPoints*(nbControlPoints-1)/2);
			// loop through the integration points
			for (l=0;l<numericalIntegrationStiffnessDataArray.size();++l) {
				Coord dpp[3],point[3];
				// the barycentric coordinate
				highOrderTriangleGeo->computeNodalValueDerivatives(i,numericalIntegrationStiffnessDataArray[l].integrationPoint, x,point);

				for(j=0; j<3; ++j){
					m=edgesInTriangleArray[j][0];
					n=edgesInTriangleArray[j][1];

					dpp[j]=point[n]-point[m];

				}
				if (decompositionMethod==QR_DECOMPOSITION)
				{
					/// perform QR decomposition
					computeQRRotation(S,dpp);
					R=S.transposed()*tetinfo->integrationPointsRestRotationArray[l];

				} else if (decompositionMethod==POLAR_DECOMPOSITION_MODIFIED) {

					S[0]=dpp[0];
					S[1]=dpp[1];
					S[2]=dpp[2];
					helper::Decompose<Real>::polarDecomposition( S, R );
					R=R.transposed()*tetinfo->integrationPointsRestRotationArray[l];
				} else {
					R.identity();
				}
				const std::vector<Vec3> & weightArray=numericalIntegrationStiffnessDataArray[l].weightArray;
				const Mat3x4  &edgeStiff1 =tetinfo->integrationPointsStiffnessVector[2*l];
				const Mat3x4  &edgeStiff2 =tetinfo->integrationPointsStiffnessVector[2*l+1];

				for (rank=0,j=0; j<nbControlPoints; ++j) {
					v0 = indexArray[j];
					for ( k=j+1; k<nbControlPoints; ++k,++rank) {

						const Vec3  & coeffVec1=weightArray[2*rank];
						const Vec3  & coeffVec2=weightArray[2*rank+1];
						Vec4 res=edgeStiff1.multTranspose(coeffVec1)+edgeStiff2.multTranspose(coeffVec2);
						Mat2x2 stiffness= Mat2x2((const Real *) &res[0]);
						//			Vec4 res=edgeStiffness[0]*coeffVec1+edgeStiffness[1]*coeffVec2;


						// loop through the integration points
						v1 = indexArray[k];
						dpos=x[v0]-x[v1];
						// displacement in the rest configuration
						dpos=R.transposed()*dpos-(x0[v0]-x0[v1]);
						// force on first vertex in the rest configuration
						force[k]-=R*stiffness*dpos;
						// force on second vertex in the rest configuration
						force[j]+=R*stiffness.multTranspose(dpos);
						if (mparams->implicit()) {
							// implicit scheme : need to store the rotated tensor
							Mat2x2 mat=R*stiffness*R.transposed();
							stiffnessArray[rank]+=mat;
						}

					}

				}
			}
//	tetinfo->rotatedStiffnessVector.clear();
//	tetinfo->rotatedStiffnessVector=stiffnessArray;
			for (j=0; j<nbControlPoints; ++j) {
				f[indexArray[j]]+=R*force[j];
			}


		} else {
			const  helper::vector<Mat2x2> &stiffnessArray=getStiffnessArray(i,tetinfo);
			assert(stiffnessArray.size()==nbControlPoints*(nbControlPoints-1)/2);
			if (decompositionMethod==LINEAR_ELASTIC) {

				for (j=0; j<nbControlPoints; ++j)
				{
					dp[j]=x[indexArray[j]]-x0[indexArray[j]];
				}
				// loop over each entry in the stiffness vector of size nbControlPoints*(nbControlPoints-1)/2
				for (rank=0,j=0; j<nbControlPoints; ++j) {
					v0 = indexArray[j];
					for ( k=j+1; k<nbControlPoints; ++k,++rank) {
						v1 = indexArray[k];
						dpos=dp[j]-dp[k];
						//		if ((i==0))
						//				std::cerr<<"stiffness["<<j<<","<<k<<"]="<<tetinfo->stiffnessVector[rank]<<" v0="<<v0<<" v1="<<v1<<std::endl;
						//					 f[v1]-=tetinfo->stiffnessVector[rank]*dp[j];
						//					 f[v0]-=tetinfo->stiffnessVector[rank].multTranspose(dp[k]);
						f[v1]-=stiffnessArray[rank]*dpos;
						f[v0]+=stiffnessArray[rank].multTranspose(dpos);
					}
				}
			} else 	{
				Mat2x2 deformationGradient,S,R;
				Coord dpp[3];
				for (j=0; j<3; ++j)
				{
					dpp[j]=x[tetinfo->v[edgesInTriangleArray[j][1]]]-x[tetinfo->v[edgesInTriangleArray[j][0]]];
				}
				if (decompositionMethod==POLAR_DECOMPOSITION)
				{
					// compute the deformation gradient
					// deformation gradient = sum of tensor product between vertex position and shape vector
					// optimize by using displacement with first vertex
					sv=tetinfo->shapeVector[1];


					for (k=0; k<2; ++k)
					{
						for (l=0; l<2; ++l)
						{
							deformationGradient[k][l]=dpp[0][k]*sv[l];
						}
					}
					for (j=1; j<2; ++j)
					{
						sv=tetinfo->shapeVector[j+1];
						for (k=0; k<2; ++k)
						{
							for (l=0; l<2; ++l)
							{
								deformationGradient[k][l]+=dpp[j][k]*sv[l];
							}
						}
					}
					// polar decomposition of the transformation
					helper::Decompose<Real>::polarDecomposition(deformationGradient,R);
				}
				else if (decompositionMethod==QR_DECOMPOSITION)
				{

					/// perform QR decomposition
					computeQRRotation(S,dpp);
					R=S.transposed()*tetinfo->restRotation;

				} else if (decompositionMethod==POLAR_DECOMPOSITION_MODIFIED) {

					S[0]=dpp[0];
					S[1]=dpp[1];
					S[2]=dpp[2];
					helper::Decompose<Real>::polarDecomposition( S, R );
					R=R.transposed()*tetinfo->restRotation;
				} 

				//	std::cerr<<"rotation= "<<R<<std::endl;
				//	R.identity();
				// store transpose of rotation
				tetinfo->rotation=R.transposed();
				std::fill(force.begin(),force.end(),Coord());
				if (mparams->implicit()) {
					tetinfo->rotatedStiffnessVector.clear();
				}
				// loop over each entry in the stiffness vector of size nbControlPoints*(nbControlPoints-1)/2
				for (rank=0,j=0; j<nbControlPoints; ++j) {
					v0 = indexArray[j];
					for ( k=j+1; k<nbControlPoints; ++k,++rank) {
						v1 = indexArray[k];
						dpos=x[v0]-x[v1];
						// displacement in the rest configuration
						dpos=tetinfo->rotation*dpos-(x0[v0]-x0[v1]);

						// force on first vertex in the rest configuration
						force[k]-=stiffnessArray[rank]*dpos;
						// force on second vertex in the rest configuration
						force[j]+=stiffnessArray[rank].multTranspose(dpos);
						if (mparams->implicit()) {
							// implicit scheme : need to store the rotated tensor
							Mat2x2 mat=R*stiffnessArray[rank]*tetinfo->rotation;
							tetinfo->rotatedStiffnessVector.push_back(mat);
						}
					}
				}
				for (j=0; j<nbControlPoints; ++j) {
					f[indexArray[j]]+=R*force[j];
				}

			}
		}

	}

    updateMatrix=true; // next time assemble the matrix
    triangleInfo.endEdit();

    dataF.endEdit();

}


template <class DataTypes>
void HighOrderTriangularCorotationalFEMForceField<DataTypes>::addDForce(const sofa::core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv&   datadF , const DataVecDeriv&   datadX )
{
    VecDeriv& df       = *(datadF.beginEdit());
    const VecCoord& dx =   datadX.getValue()  ;
    Real kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());
    const helper::vector<TriangleRestInformation>& triangleInf = triangleInfo.getValue();
	Coord dpos;
    size_t i,j,k,v0,v1,rank;
    size_t nbTriangles=_topology->getNbTriangles();
    
	HighOrderDegreeType degree=highOrderTriangleGeo->getTopologyContainer()->getDegree();
	size_t nbControlPoints=(degree+1)*(degree+2)/2;
	
    const TriangleRestInformation *tetinfo;
	

	for(i=0; i<nbTriangles; i++ )
	{
		tetinfo=&triangleInf[i];
		const HighOrderTriangleSetTopologyContainer::VecPointID &indexArray=
			highOrderTriangleGeo->getTopologyContainer()->getGlobalIndexArrayOfControlPoints(i);
		const  helper::vector<Mat2x2> &stiffnessArray=getRotatedStiffnessArray(i,tetinfo);
		// create a local buffer to limit access to the df array
		sofa::helper::vector<Deriv> dforce;

		nbControlPoints=indexArray.size();
		dforce.resize(nbControlPoints);
		assert(stiffnessArray.size()==nbControlPoints*(nbControlPoints-1)/2);
		// loop over each entry in the stiffness vector of size nbControlPoints*(nbControlPoints+1)/2
		for (rank=0,j=0; j<nbControlPoints; ++j) {
			v0 = indexArray[j];
			for ( k=j+1; k<nbControlPoints; ++k,++rank) {
				v1 = indexArray[k];
				dpos=dx[v0]-dx[v1];
				dforce[k]-=stiffnessArray[rank]*dpos*kFactor;
				dforce[j]+=stiffnessArray[rank].multTranspose(dpos*kFactor);
			}
		}
		for (j=0; j<nbControlPoints; ++j) {
			df[indexArray[j]]+=dforce[j];
		}
	
	}

    datadF.endEdit();
}

template<class DataTypes>
SReal HighOrderTriangularCorotationalFEMForceField<DataTypes>::getPotentialEnergy(const core::MechanicalParams*, const DataVecCoord&) const
{
    serr << "ERROR("<<this->getClassName()<<"): getPotentialEnergy( const MechanicalParams*, const DataVecCoord& ) not implemented." << sendl;
    return 0.0;
}


template<class DataTypes>
void HighOrderTriangularCorotationalFEMForceField<DataTypes>::updateLameCoefficients()
{
    lambda= d_youngModulus.getValue()*d_poissonRatio.getValue()/(1-d_poissonRatio.getValue()*d_poissonRatio.getValue());
    mu = d_youngModulus.getValue()/((1+d_poissonRatio.getValue()));

}


template<class DataTypes>
void HighOrderTriangularCorotationalFEMForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
#ifndef SOFA_NO_OPENGL
    if (!vparams->displayFlags().getShowForceFields()) return;
    if (!this->mstate) return;

    if (vparams->displayFlags().getShowWireFrame())
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);


    if (vparams->displayFlags().getShowWireFrame())
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
#endif /* SOFA_NO_OPENGL */
}

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_HighOrderTriangularCorotationalFEMForceField_INL
