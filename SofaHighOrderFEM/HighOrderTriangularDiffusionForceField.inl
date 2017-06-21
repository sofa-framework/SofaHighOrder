
#ifndef SOFA_COMPONENT_FORCEFIELD_HIGHORDERTRIANGULARCOROTATIONALFEMFORCEFIELD_INL
#define SOFA_COMPONENT_FORCEFIELD_HIGHORDERTRIANGULARCOROTATIONALFEMFORCEFIELD_INL

#include "HighOrderTriangularDiffusionForceField.h"
#include <sofa/core/visual/VisualParams.h>
#include <fstream> // for reading the file
#include <iostream> //for debugging
#include <sofa/helper/gl/template.h>
#include <SofaBaseTopology/TopologyData.inl>
#include <HighOrderTriangleSetGeometryAlgorithms.h>
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

typedef EdgesInTriangle		EdgesInTriangle;

const unsigned int edgesInTriangleArray[3][2] = {{0,1}, {0,2}, {1,2}};


template< class DataTypes>
void HighOrderTriangularDiffusionForceField<DataTypes>::weightArrayPointer::allocate(size_t degree) {
	if (degree==2) {
		weightArrayQuadratic=boost::make_shared<Mat15x3>();
		

	} else 	if (degree==3) {
		weightArrayCubic=boost::make_shared<Mat45x3>();
		
	} else 	if (degree==4) {
		weightArrayQuartic=boost::make_shared<Mat105x3>();
		
	} else 	if (degree==5) {
		weightArrayQuintic=boost::make_shared<Mat210x3>();
		
	}
}


template< class DataTypes>
void HighOrderTriangularDiffusionForceField<DataTypes>::FTCFTriangleHandler::applyCreateFunction(unsigned int triangleIndex,
        TriangleRestInformation &my_tinfo,
        const Triangle &,
        const sofa::helper::vector<unsigned int> &,
        const sofa::helper::vector<double> &)
{
	if (ff)
	{
		const helper::vector< Triangle > &triangleArray=ff->_topology->getTriangles() ;
		HighOrderTriangleSetTopologyContainer *container=ff->highOrderTrianGeo->getTopologyContainer();
		HighOrderDegreeType degree=container->getDegree();
		size_t nbControlPoints=(degree+1)*(degree+2)/2;
		size_t nbStiffnessEntries=nbControlPoints*(nbControlPoints-1)/2;
		if (my_tinfo.stiffnessVector.size()!=nbStiffnessEntries) {
			my_tinfo.stiffnessVector.resize(nbStiffnessEntries);
		}
		// set array to zero
		std::fill(my_tinfo.stiffnessVector.begin(),my_tinfo.stiffnessVector.end(),Real());

		//		const std::vector< Edge> &edgeArray=ff->_topology->getEdges() ;
		size_t i,j,k,l,m,n;

		Vec2 point[3];


		Vec3 edgeStiffness;  // the off-diagonal 3x3 block matrices that makes the 12x12 linear elastic matrix

		 const typename HighOrderTriangularDiffusionForceField<DataTypes>::MechanicalTypes::VecCoord  &	restPosition=ff->mechanicalObject->read(core::ConstVecCoordId::restPosition())->getValue();
		// now computed the stiffness for the HighOrder Triangle
		sofa::helper::vector<TriangleIndexVector> tbiArray;

		tbiArray=ff->highOrderTrianGeo->getTopologyContainer()->getTriangleIndexArray();

		size_t rank=0;

		///describe the indices of the 3 triangle vertices
		const Triangle &t= triangleArray[triangleIndex];
		//    BaseMeshTopology::EdgesInTriangle te=ff->_topology->getEdgesInTriangle(triangleIndex);


		// store the point position
		for(j=0; j<3; ++j)
			point[j]=(restPosition)[t[j]];

		if ((ff->integrationMethod==HighOrderTriangularDiffusionForceField<DataTypes>::AFFINE_ELEMENT_INTEGRATION) || 
			((ff->d_forceAffineAssemblyForAffineElements.getValue()) && ((ff->highOrderTrianGeo->isBezierTriangleAffine(triangleIndex,restPosition  ))))) {


				helper::system::thread::ctime_t startUpdateMat=helper::system::thread::CTime::getTime();
				ff->computeTriangleStiffnessEdgeMatrix(point,edgeStiffness);
				if (degree==1) {
				for (rank=0;rank<nbStiffnessEntries;rank++) {
					my_tinfo.stiffnessVector[rank]+=  edgeStiffness[rank];
				}
			} else if (degree==2) {
				Vec15 res=(*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayQuadratic))*edgeStiffness;

				for (rank=0;rank<nbStiffnessEntries;rank++) {
					my_tinfo.stiffnessVector[rank]+=  res[rank];
				}
			} else if (degree==3) {
				Vec45 res=(*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayCubic))*edgeStiffness;

				for (rank=0;rank<nbStiffnessEntries;rank++) {
					my_tinfo.stiffnessVector[rank]+=  res[rank];
				}
			} else if (degree==4) {
				Vec105 res=(*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayQuartic))*edgeStiffness;

				for (rank=0;rank<nbStiffnessEntries;rank++) {
					my_tinfo.stiffnessVector[rank]+=  res[rank];
				}
			} else if (degree==5) {
				Vec210 res=(*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayQuintic))*edgeStiffness;

				for (rank=0;rank<nbStiffnessEntries;rank++) {
					my_tinfo.stiffnessVector[rank]+=  res[rank];
				}

			} else {
				for (rank=0;rank<nbStiffnessEntries;rank++) {
					my_tinfo.stiffnessVector[rank]+= dot(edgeStiffness,ff->affineStiffnessCoefficientArray[rank]);
				}
			}


			if (ff->f_printLog.getValue()) {
				helper::system::thread::ctime_t endUpdateMat=helper::system::thread::CTime::getTime();
				ff->totalUpdateMat+=endUpdateMat-startUpdateMat;
			}
		} else if (ff->integrationMethod==HighOrderTriangularDiffusionForceField<DataTypes>::STANDARD_INTEGRATION){
			sofa::defaulttype::Vec<3,Real> bc;
			helper::system::thread::ctime_t startUpdateMat=helper::system::thread::CTime::getTime();
			// loop through the integration points
			for (i=0;i<ff->numericalIntegrationStiffnessDataArray.size();++i) {

				// the barycentric coordinate
				bc=ff->numericalIntegrationStiffnessDataArray[i].integrationPoint;
				// Compute the local triangle by storing in point the 3 derivatives of the shape functions  
				ff->highOrderTrianGeo->computeNodalValueDerivatives(triangleIndex,bc, restPosition,point);

				Mat2x2 Jacobian,inverseJacobian;
				for (j=0;j<2;++j) {
					for (k=0;k<2;++k) {
						Jacobian[j][k]=point[j][k]-point[2][k];
					}
				}
				invertMatrix(inverseJacobian,Jacobian);
				Real jac=fabs(determinant(Jacobian))*ff->numericalIntegrationStiffnessDataArray[i].integrationWeight;

				helper::vector<Vec2> SVArray;
				for (j=0;j<nbControlPoints;j++) {
					Vec2 sv=inverseJacobian*ff->numericalIntegrationStiffnessDataArray[i].coefficientArray[j];
					SVArray.push_back(sv);
				}
				Real coeff=ff->d_diffusivity.getValue()*jac;
				for (rank=0,j=0;j<nbControlPoints;j++) {

					for (k=j+1;k<nbControlPoints;k++,rank++) {
						if (ff->diffusionSymmetry==ISOTROPIC) {
							my_tinfo.stiffnessVector[rank]+=dot(SVArray[j], SVArray[k])*coeff;
						} else {
							my_tinfo.stiffnessVector[rank]+=dot(SVArray[j], ff->diffusionTensor*SVArray[k])*jac;
						}
					}
				}
			}
			if (ff->f_printLog.getValue()) {
				helper::system::thread::ctime_t endUpdateMat=helper::system::thread::CTime::getTime();
				ff->totalUpdateMat+=endUpdateMat-startUpdateMat;
			}

		} else if (ff->integrationMethod==HighOrderTriangularDiffusionForceField<DataTypes>::NUMERICAL_INTEGRATION_2){
			helper::system::thread::ctime_t startUpdateMat=helper::system::thread::CTime::getTime();
			sofa::defaulttype::Vec<3,Real> bc;

			// loop through the integration points
			for (i=0;i<ff->numericalIntegrationStiffnessDataArray.size();++i) {

				// the barycentric coordinate
				bc=ff->numericalIntegrationStiffnessDataArray[i].integrationPoint;
				// Compute the local triangle by storing in point the 3 derivatives of the shape functions  
				ff->highOrderTrianGeo->computeNodalValueDerivatives(triangleIndex,bc, restPosition,point);
				// compute the edge stiffness associated with that local triangle
				ff->computeTriangleStiffnessEdgeMatrix(point,edgeStiffness);

				// compute the stiffness matrix for all pairs of control points 
				for (rank=0,j=0;j<nbControlPoints;j++) {

					for (k=j+1;k<nbControlPoints;k++,rank++) {

						const Mat3x3  & coeffMatrix=ff->numericalIntegrationStiffnessDataArray[i].weightArray[rank];
						/// add edge stiffness
						for(l=0;l<3;++l) {
							m=edgesInTriangleArray[l][0];
							n=edgesInTriangleArray[l][1];
							if (coeffMatrix[m][n]!=0) {
								my_tinfo.stiffnessVector[rank]+= edgeStiffness[l]*coeffMatrix[m][n];
							}	
						}
					}
				}
			}
			if (ff->f_printLog.getValue()) {
				helper::system::thread::ctime_t endUpdateMat=helper::system::thread::CTime::getTime();
				ff->totalUpdateMat+=endUpdateMat-startUpdateMat;
			}
		} else if (ff->integrationMethod==HighOrderTriangularDiffusionForceField<DataTypes>::NUMERICAL_INTEGRATION){
			helper::system::thread::ctime_t startUpdateMat=helper::system::thread::CTime::getTime();
			sofa::defaulttype::Vec<3,Real> bc;

			// loop through the integration points
			for (i=0;i<ff->numericalIntegrationStiffnessDataArray.size();++i) {

				// the barycentric coordinate
				bc=ff->numericalIntegrationStiffnessDataArray[i].integrationPoint;
				// Compute the local triangle by storing in point the 3 derivatives of the shape functions  
				ff->highOrderTrianGeo->computeNodalValueDerivatives(triangleIndex,bc, restPosition,point);
				// compute the edge stiffness associated with that local triangle
				ff->computeTriangleStiffnessEdgeMatrix(point,edgeStiffness);
				if (degree==2) {
					Vec15 res=(*(ff->numericalIntegrationStiffnessDataArray[i].arrayPointer.weightArrayQuadratic))*edgeStiffness;
					for (rank=0;rank<nbStiffnessEntries;rank++) {
						my_tinfo.stiffnessVector[rank]+=  res[rank];
					}
				} 	else if (degree==3) {
					Vec45 res=(*(ff->numericalIntegrationStiffnessDataArray[i].arrayPointer.weightArrayCubic))*edgeStiffness;
					for (rank=0;rank<nbStiffnessEntries;rank++) {
						my_tinfo.stiffnessVector[rank]+=  res[rank];
					}
				} 	else if (degree==4) {
					Vec105 res=(*(ff->numericalIntegrationStiffnessDataArray[i].arrayPointer.weightArrayQuartic))*edgeStiffness;
					for (rank=0;rank<nbStiffnessEntries;rank++) {
						my_tinfo.stiffnessVector[rank]+=  res[rank];
					}
				} 	else if (degree==5) {
					Vec210 res=(*(ff->numericalIntegrationStiffnessDataArray[i].arrayPointer.weightArrayQuintic))*edgeStiffness;
					for (rank=0;rank<nbStiffnessEntries;rank++) {
						my_tinfo.stiffnessVector[rank]+=  res[rank];
					}
				} else {
					// compute the stiffness matrix for all pairs of control points 
					for (rank=0;rank<nbStiffnessEntries;rank++) {
						const Vec3  & coeffVec=ff->numericalIntegrationStiffnessDataArray[i].weightVectorizedArray[rank];
						my_tinfo.stiffnessVector[rank]+= dot(edgeStiffness,coeffVec);
					}
				}		
			}
			if (ff->f_printLog.getValue()) {
				helper::system::thread::ctime_t endUpdateMat=helper::system::thread::CTime::getTime();
				ff->totalUpdateMat+=endUpdateMat-startUpdateMat;
			}
		}
				
#ifdef _DEBUG
		if (ff->f_printLog.getValue()) {
			std::cerr<< " diffusivity="<<ff->d_diffusivity.getValue() << std::endl;


			for (rank=0,l=0;l<tbiArray.size();++l)
			{
				for (m=l+1;m<tbiArray.size();++m,++rank)
				{
					std::cerr<< "Stiffness entry ["<<(unsigned int)tbiArray[l][0]<<" "<< (unsigned int)tbiArray[l][1]<<" "<< (unsigned int)tbiArray[l][2]<<" ]["<<
						(unsigned int)tbiArray[m][0]<<" "<< (unsigned int)tbiArray[m][1]<<" "<< (unsigned int)tbiArray[m][2]<<" "<<"]="<<my_tinfo.stiffnessVector[rank]<<std::endl;

				}
			}
		}
#endif
	}
}



template <class DataTypes> HighOrderTriangularDiffusionForceField<DataTypes>::HighOrderTriangularDiffusionForceField()
    : triangleInfo(initData(&triangleInfo, "triangleInfo", "Internal triangle data"))
    , updateMatrix(true)
    , d_diffusivity(initData(&d_diffusivity,(Real)1000.,"diffusivity","diffusivity for isotropic diffusion"))
	, d_anisotropy(initData(&d_anisotropy,std::string("isotropy"),"anisotropy","\"isotropy\" or \"transverseIsotropy\" or \"orthotropy\" as an anisotropy of diffusion"))
	, numericalIntegrationOrder( initData(&numericalIntegrationOrder,(size_t)2,"integrationOrder","The order of integration for numerical integration"))
	, d_integrationMethod( initData(&d_integrationMethod,std::string("analytical"),"integrationMethod","\"analytical\" if closed form expression for affine element, \"numerical\" if numerical integration is chosen,  \"standard\" if standard integration is chosen"))
	, d_anisotropyParameter( initData(&d_anisotropyParameter,ParameterArray(),"anisotropyParameters","diffusivity parameters if the diffusion is anisotropic "))
	, d_anisotropyDirection( initData(&d_anisotropyDirection,AnisotropyDirectionArray(),"anisotropyDirection","Directions of anisotropy"))
	, numericalIntegrationMethod( initData(&numericalIntegrationMethod,(size_t)0,"numericalIntegrationMethod","The type of numerical integration method chosen"))
    , d_assemblyTime(initData(&d_assemblyTime,(Real)0,"assemblyTime","the time spent in assembling the stiffness matrix. Only updated if printLog is set to true"))
	 , d_forceAffineAssemblyForAffineElements(initData(&d_forceAffineAssemblyForAffineElements,true,"forceAffineAssemblyForAffineElements","if true affine triangles are always assembled with the closed form formula, Otherwise use the method defined in integrationMethod"))
    , triangleHandler(NULL)
{
    triangleHandler = new FTCFTriangleHandler(this,&triangleInfo);
}

template <class DataTypes> HighOrderTriangularDiffusionForceField<DataTypes>::~HighOrderTriangularDiffusionForceField()
{
    if (triangleHandler) delete triangleHandler;
}


template <class DataTypes> void HighOrderTriangularDiffusionForceField<DataTypes>::init()
{
    //	serr << "initializing HighOrderTriangularDiffusionForceField" << sendl;
    this->Inherited::init();

    _topology = this->getContext()->getMeshTopology();
	this->getContext()->get(highOrderTrianGeo);

    if ((_topology->getNbTriangles()==0) || (!highOrderTrianGeo))
    {
        serr << "ERROR(HighOrderTriangularDiffusionForceField): object must have a Triangular Set Topology and a HighOrderTriangleSetGeometryAlgorithms component "<<sendl;
        return;
    }
   
	if (d_integrationMethod.getValue() == "analytical")
        integrationMethod= AFFINE_ELEMENT_INTEGRATION;
    else if (d_integrationMethod.getValue() == "numerical") 
		integrationMethod= NUMERICAL_INTEGRATION;
	else if (d_integrationMethod.getValue() == "numerical2") 
		integrationMethod= NUMERICAL_INTEGRATION_2;
	else if (d_integrationMethod.getValue() == "standard") 
		integrationMethod= STANDARD_INTEGRATION;
    else
    {
        serr << "cannot recognize method "<< d_integrationMethod.getValue() << ". Must be either \"analytical\" or \"numerical\"  or \"standard\"" << sendl;
    }
	if (d_anisotropy.getValue() == "isotropy")
        diffusionSymmetry= ISOTROPIC;
    else if (d_anisotropy.getValue() == "transverseIsotropy") 
        diffusionSymmetry= TRANSVERSE_ISOTROPIC;
    else if (d_anisotropy.getValue() == "orthotropy") 
        diffusionSymmetry= ORTHOTROPIC;
    else
    {
        serr << "cannot recognize anisotropy type "<< d_anisotropy.getValue() << ". Must be either \"isotropy\" or \"transverseIsotropy\"  or \"orthotropy\"" << sendl;
    }

    helper::vector<TriangleRestInformation>& triangleInf = *(triangleInfo.beginEdit());
    triangleInf.resize(_topology->getNbTriangles());


	/// Get the mechanical object containing the mesh position in 3D
	sofa::core::objectmodel::Tag mechanicalTag(m_tagMeshMechanics.getValue());
	this->getContext()->get(mechanicalObject, mechanicalTag,sofa::core::objectmodel::BaseContext::SearchUp);
	if (mechanicalObject==NULL)
	{
		serr<<"ERROR(HighOrderTriangularDiffusionForceField): cannot find the mechanical object."<<sendl;
		std::cout<<"mechanicalObj = "<<mechanicalObject<<std::endl;
		return;
	}

	helper::system::thread::ctime_t startAssembly=helper::system::thread::CTime::getTime();
	totalUpdateMat=0;

	size_t i;

	topology::HighOrderDegreeType degree=highOrderTrianGeo->getTopologyContainer()->getDegree();
	size_t nbControlPoints=(degree+1)*(degree+2)/2;
	sofa::helper::vector<topology::TriangleIndexVector> tbiArray;
	tbiArray=highOrderTrianGeo->getTopologyContainer()->getTriangleIndexArray();
	topology::TriangleIndexVector tbi1,tbi2;

	if (degree==1) 
		integrationMethod= AFFINE_ELEMENT_INTEGRATION;

	computeDiffusivityTensor();
	if (degree>1) {
		if ((integrationMethod== AFFINE_ELEMENT_INTEGRATION) || (d_forceAffineAssemblyForAffineElements.getValue()))
		{
			// precompute the coefficients for handling affine elements
			if (degree<6) {
				affineStiffnessCoefficientPreStoredArray.allocate(degree);
			}
			affineStiffnessCoefficientArray.clear();
			std::vector<Real> coeffArray(3);
			Mat3x3 coeffMatrix;
			size_t j,k,l,m,n,rank;
			for (rank=0,j=0;j<nbControlPoints;j++) {
				tbi1=tbiArray[j];
				for (k=j+1;k<nbControlPoints;k++,rank++) {
					tbi2=tbiArray[k];
					coeffMatrix=highOrderTrianGeo->getAffineStiffnessCoefficientMatrix(tbi1,tbi2);
					// substract the diagonal terms such that only edge stiffness are used
					for(l=0; l<3; ++l){
						for(m=l+1; m<3; ++m){
							coeffMatrix[l][m]+= coeffMatrix[m][l]-(coeffMatrix[l][l]+coeffMatrix[m][m]);						

						}
					}
					if (degree==2) {
						for(l=0; l<3; ++l){
							m=edgesInTriangleArray[l][0];
							n=edgesInTriangleArray[l][1];

							(*(affineStiffnessCoefficientPreStoredArray.weightArrayQuadratic))[rank][l]=coeffMatrix[m][n];

						}
					} else if (degree==3) {
						for(l=0; l<3; ++l){
							m=edgesInTriangleArray[l][0];
							n=edgesInTriangleArray[l][1];

							(*(affineStiffnessCoefficientPreStoredArray.weightArrayCubic))[rank][l]=coeffMatrix[m][n];

						}
					} else if (degree==4) {
						for(l=0; l<3; ++l){
							m=edgesInTriangleArray[l][0];
							n=edgesInTriangleArray[l][1];

							(*(affineStiffnessCoefficientPreStoredArray.weightArrayQuartic))[rank][l]=coeffMatrix[m][n];

						}
					} else if (degree==5) {
						for(l=0; l<3; ++l){
							m=edgesInTriangleArray[l][0];
							n=edgesInTriangleArray[l][1];

							(*(affineStiffnessCoefficientPreStoredArray.weightArrayQuintic))[rank][l]=coeffMatrix[m][n];

						}

					} else {
						Vec3 coeff;
						for(l=0; l<3; ++l){
							m=edgesInTriangleArray[l][0];
							n=edgesInTriangleArray[l][1];
							coeff[l]=coeffMatrix[m][n];
						}
						affineStiffnessCoefficientArray.push_back(coeff);
					}
				}
			}
		}
		if ((integrationMethod== NUMERICAL_INTEGRATION) || (integrationMethod== NUMERICAL_INTEGRATION_2))
		{
			numericalIntegrationStiffnessDataArray.clear();
			/// get value of integration points0
			topology::NumericalIntegrationDescriptor<Real,3> &nid=highOrderTrianGeo->getTriangleNumericalIntegrationDescriptor();
			typename topology::NumericalIntegrationDescriptor<Real,3>::QuadraturePointArray qpa=nid.getQuadratureMethod((typename topology::NumericalIntegrationDescriptor<Real,3>::QuadratureMethod)numericalIntegrationMethod.getValue(),
				numericalIntegrationOrder.getValue());
			size_t i,j,k,l,m;
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
				size_t n;

				std::vector<Vec3> shapeFunctionDerivativeArray;
				for(j=0;j<tbiArray.size();++j) {
					Vec3 deriv=highOrderTrianGeo->computeShapeFunctionDerivatives(tbiArray[j],qp.first);
					shapeFunctionDerivativeArray.push_back(deriv);
					Vec2 der(deriv[0]-deriv[2],deriv[1]-deriv[2]);
					nimd.coefficientArray.push_back(der);
				}
				size_t rank;
				for(rank=0,j=0;j<tbiArray.size();++j) {
					for(k=j+1;k<tbiArray.size();++k,++rank) {
						coeffMatrix=dyad(shapeFunctionDerivativeArray[j],shapeFunctionDerivativeArray[k])*2*weight;
						for(l=0; l<3; ++l){
							for(m=l+1; m<3; ++m){
								coeffMatrix[l][m]+= coeffMatrix[m][l]-(coeffMatrix[l][l]+coeffMatrix[m][m]);
							}
						}
						if (integrationMethod== NUMERICAL_INTEGRATION) { 
							if (degree>5)  {
								Vec3 coeffVec;
								for(l=0; l<3; ++l){
									m=edgesInTriangleArray[l][0];
									n=edgesInTriangleArray[l][1];
									coeffVec[l]=coeffMatrix[m][n];								
								}
								nimd.weightVectorizedArray.push_back(coeffVec);

							} else {
								if (degree==2) {

									for(l=0; l<3; ++l){
										m=edgesInTriangleArray[l][0];
										n=edgesInTriangleArray[l][1];
										(*(nimd.arrayPointer.weightArrayQuadratic))[rank][l]=coeffMatrix[m][n];
									}
								} else 	if (degree==3) {

									for(l=0; l<3; ++l){
										m=edgesInTriangleArray[l][0];
										n=edgesInTriangleArray[l][1];
										(*(nimd.arrayPointer.weightArrayCubic))[rank][l]=coeffMatrix[m][n];
									}
								} else 	if (degree==4) {


									for(l=0; l<3; ++l){
										m=edgesInTriangleArray[l][0];
										n=edgesInTriangleArray[l][1];
										(*(nimd.arrayPointer.weightArrayQuartic))[rank][l]=coeffMatrix[m][n];
									}
								} else 	if (degree==5) {

									for(l=0; l<3; ++l){
										m=edgesInTriangleArray[l][0];
										n=edgesInTriangleArray[l][1];
										(*(nimd.arrayPointer.weightArrayQuintic))[rank][l]=coeffMatrix[m][n];
									}
								}
							}	
						} else 	if (integrationMethod== NUMERICAL_INTEGRATION_2) { 
							nimd.weightArray.push_back(coeffMatrix);
						}

					}
				}

				numericalIntegrationStiffnessDataArray.push_back(nimd);
			}
		}
		if (integrationMethod== STANDARD_INTEGRATION) 
		{
			numericalIntegrationStiffnessDataArray.clear();
			/// get value of integration points0
			topology::NumericalIntegrationDescriptor<Real,3> &nid=highOrderTrianGeo->getTriangleNumericalIntegrationDescriptor();
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
					Vec3 deriv=highOrderTrianGeo->computeShapeFunctionDerivatives(tbiArray[j],qp.first);
					Vec2 der(deriv[0]-deriv[2],deriv[1]-deriv[2]);
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
void HighOrderTriangularDiffusionForceField<DataTypes>::updateTopologyInformation()
{
    int i;
    unsigned int j;

    int nbTriangles=_topology->getNbTriangles();

    TriangleRestInformation *tetinfo;

    helper::vector<typename HighOrderTriangularDiffusionForceField<DataTypes>::TriangleRestInformation>& triangleInf = *(triangleInfo.beginEdit());
   


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
void HighOrderTriangularDiffusionForceField<DataTypes>::computeDiffusivityTensor() 											
{
	const helper::vector<Real> & anisotropyParameter=d_anisotropyParameter.getValue();
	const helper::vector<Vec2> & anisotropyDirection=d_anisotropyDirection.getValue();
	if (diffusionSymmetry==ISOTROPIC) {

	} else if (diffusionSymmetry==TRANSVERSE_ISOTROPIC) {
		assert(anisotropyParameter.size()>=1);
		assert(anisotropyDirection.size()>=1);
		diffusionTensor.identity();

		diffusionTensor+=dyad(anisotropyDirection[0],anisotropyDirection[0])*(anisotropyParameter[0]-1);
		diffusionTensor*=d_diffusivity.getValue();
	}
}

template<class DataTypes>
void HighOrderTriangularDiffusionForceField<DataTypes>::computeTriangleStiffnessEdgeMatrix(const Vec2 point[3],
																							   Vec3 &edgeStiffness)
{
	Vec2 shapeVector[3];
	/// compute 6 times the rest volume
	Real volume=areaProduct(point[1]-point[0],point[2]-point[0]);
	/// store the rest volume
//	my_tinfo.restVolume=volume/6;

	size_t j,k,l;
	// store shape vectors at the rest configuration
	for(j=0; j<3; ++j)
	{
		shapeVector[j]=ortho(point[(j+2)%3]-point[(j+1)%3])/volume;		
	}
	if (diffusionSymmetry==ISOTROPIC) {
		Real diff=d_diffusivity.getValue()*fabs(volume)/2;


		/// compute the edge stiffness of the linear elastic material
		for(j=0; j<3; ++j)
		{
			k=edgesInTriangleArray[j][0];
			l=edgesInTriangleArray[j][1];
			// the linear stiffness matrix using shape vectors and diffusivity coefficients
			edgeStiffness[j]= -diff*dot(shapeVector[l],shapeVector[k]);
		}
	} else {
		Mat2x2 diff=diffusionTensor*fabs(volume)/2;
		for(j=0; j<3; ++j)
		{
			k=edgesInTriangleArray[j][0];
			l=edgesInTriangleArray[j][1];
			edgeStiffness[j]= -dot(shapeVector[l],diff*shapeVector[k]);
		}
	}
}

template <class DataTypes>
const helper::vector<typename HighOrderTriangularDiffusionForceField<DataTypes>::Real> & 
	HighOrderTriangularDiffusionForceField<DataTypes>::getStiffnessArray(
	const size_t i,
	const typename HighOrderTriangularDiffusionForceField<DataTypes>::TriangleRestInformation *restTrian)
{
	return(restTrian->stiffnessVector);
}

template <class DataTypes>
void HighOrderTriangularDiffusionForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* mparams /* PARAMS FIRST */, 
																	 DataVecDeriv &  dataF, const DataVecCoord &  dataX , const DataVecDeriv & /*dataV*/ )
{

    VecDeriv& f        = *(dataF.beginEdit());
    const VecCoord& x  =   dataX.getValue()  ;
	
    size_t i,j,k,l,v0,v1,rank,p;
	size_t nbTriangles=_topology->getNbTriangles();
	HighOrderDegreeType degree=highOrderTrianGeo->getTopologyContainer()->getDegree();
	size_t nbControlPoints=(degree+1)*(degree+2)/2;

	HighOrderTriangleSetTopologyContainer::VecPointID indexArray;


    if (updateTopologyInfo)
    {
        updateTopologyInformation();
    }
    helper::vector<TriangleRestInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleRestInformation *tetinfo;

	
	
	Coord dpos,sv;
		
	for(i=0; i<nbTriangles; i++ )
	{
		tetinfo=&triangleInf[i];
		const HighOrderTriangleSetTopologyContainer::VecPointID &indexArray=
			highOrderTrianGeo->getTopologyContainer()->getGlobalIndexArrayOfControlPoints(i);
		// reset force vector
		sofa::helper::vector<Coord> force;
		force.resize(nbControlPoints);
		const  helper::vector<Real> &stiffnessArray=getStiffnessArray(i,tetinfo);

		// loop over each entry in the stiffness vector of size nbControlPoints*(nbControlPoints-1)/2
		for (rank=0,j=0; j<nbControlPoints; ++j) {
			v0 = indexArray[j];
			for ( k=j+1; k<nbControlPoints; ++k,++rank) {
				v1 = indexArray[k];
				dpos=stiffnessArray[rank]*(x[v0]-x[v1]);

				force[k]-= dpos;
				force[j]+= dpos;
			}
		}
		for (j=0; j<nbControlPoints; ++j) {
			f[indexArray[j]]+=force[j];
		}
			
	}
    updateMatrix=true; // next time assemble the matrix
    triangleInfo.endEdit();

    dataF.endEdit();

}


template <class DataTypes>
void HighOrderTriangularDiffusionForceField<DataTypes>::addDForce(const sofa::core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv&   datadF , const DataVecDeriv&   datadX )
{
    VecDeriv& df       = *(datadF.beginEdit());
    const VecCoord& dx =   datadX.getValue()  ;
    Real kFactor = (Real)mparams->kFactor();
    const helper::vector<TriangleRestInformation>& triangleInf = triangleInfo.getValue();
	Coord dpos;
    size_t i,j,k,v0,v1,rank;
    size_t nbTriangles=_topology->getNbTriangles();
    
	HighOrderDegreeType degree=highOrderTrianGeo->getTopologyContainer()->getDegree();
	size_t nbControlPoints=(degree+1)*(degree+2)/2;
	
    const TriangleRestInformation *tetinfo;
	

	for(i=0; i<nbTriangles; i++ )
	{
		tetinfo=&triangleInf[i];
		const HighOrderTriangleSetTopologyContainer::VecPointID &indexArray=
			highOrderTrianGeo->getTopologyContainer()->getGlobalIndexArrayOfControlPoints(i);
		const  helper::vector<Real> &stiffnessArray=getStiffnessArray(i,tetinfo);
		// create a local buffer to limit access to the df array
		sofa::helper::vector<Deriv> dforce;

		nbControlPoints=indexArray.size();
		dforce.resize(nbControlPoints);
		
		// loop over each entry in the stiffness vector of size nbControlPoints*(nbControlPoints+1)/2
		for (rank=0,j=0; j<nbControlPoints; ++j) {
			v0 = indexArray[j];
			for ( k=j+1; k<nbControlPoints; ++k,++rank) {
				v1 = indexArray[k];
				dpos=stiffnessArray[rank]*(dx[v0]-dx[v1]);
				dforce[k]-=dpos*kFactor;
				dforce[j]+=dpos*kFactor;
			}
		}
		for (j=0; j<nbControlPoints; ++j) {
			df[indexArray[j]]+=dforce[j];
		}
	}

    datadF.endEdit();
}

template<class DataTypes>
SReal HighOrderTriangularDiffusionForceField<DataTypes>::getPotentialEnergy(const core::MechanicalParams*, const DataVecCoord&) const
{
    serr << "ERROR("<<this->getClassName()<<"): getPotentialEnergy( const MechanicalParams*, const DataVecCoord& ) not implemented." << sendl;
    return 0.0;
}




template<class DataTypes>
void HighOrderTriangularDiffusionForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
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

#endif // SOFA_COMPONENT_FORCEFIELD_HighOrderTriangularDiffusionForceField_INL
