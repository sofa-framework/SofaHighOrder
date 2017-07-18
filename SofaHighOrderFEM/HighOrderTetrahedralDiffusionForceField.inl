
#ifndef SOFA_COMPONENT_FORCEFIELD_HIGHORDERTETRAHEDRALCOROTATIONALFEMFORCEFIELD_INL
#define SOFA_COMPONENT_FORCEFIELD_HIGHORDERTETRAHEDRALCOROTATIONALFEMFORCEFIELD_INL

#include "HighOrderTetrahedralDiffusionForceField.h"
#include <sofa/core/visual/VisualParams.h>
#include <fstream> // for reading the file
#include <iostream> //for debugging
#include <sofa/helper/gl/template.h>
#include <SofaBaseTopology/TopologyData.inl>
#include <HighOrderTetrahedronSetGeometryAlgorithms.h>
#include <BezierTetrahedronSetGeometryAlgorithms.h>
#include <sofa/core/behavior/ForceField.inl>

#include <sofa/helper/decompose.h>
#include <boost/make_shared.hpp>
#include "GenericMatrixManipulator.inl"

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

typedef BaseMeshTopology::Tetra				Tetra;
typedef BaseMeshTopology::EdgesInTetrahedron		EdgesInTetrahedron;

typedef Tetra			        Tetrahedron;
typedef EdgesInTetrahedron		EdgesInTetrahedron;

const unsigned int edgesInTetrahedronArray[6][2] = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}};


template< class DataTypes>
void HighOrderTetrahedralDiffusionForceField<DataTypes>::weightArrayPointer::allocate(size_t degree) {
	if (degree==2) {
		weightArrayQuadratic=boost::make_shared<Mat45x6>();
	} else 	if (degree==3) {
		weightArrayCubic=boost::make_shared<Mat190x6>();
	} else 	if (degree==4) {
		weightArrayQuartic=boost::make_shared<Mat595x6>();
	} else 	if (degree==5) {
		weightArrayQuintic=boost::make_shared<Mat1540x6>();	
	}
}
template< class DataTypes>
void HighOrderTetrahedralDiffusionForceField<DataTypes>::weightVectorPointer::allocate(size_t degree) {
	if (degree==2) {
		resultQuadratic=boost::make_shared<Vec45>();

	} else 	if (degree==3) {
		resultCubic=boost::make_shared<Vec190>();
	} else 	if (degree==4) {
		resultQuartic=boost::make_shared<Vec595>();
	} else 	if (degree==5) {
		resultQuintic=boost::make_shared<Vec1540>();		
	}
}
template< class DataTypes>
void HighOrderTetrahedralDiffusionForceField<DataTypes>::weightVectorPointer::resetResult(size_t degree) {
	if (degree==2) {
		(*resultQuadratic).clear();
	} else 	if (degree==3) {
		(*resultCubic).clear();
	} else 	if (degree==4) {
		(*resultQuartic).clear();
	} else 	if (degree==5) {
		(*resultQuintic).clear();	
	}
}
template< class DataTypes>
void HighOrderTetrahedralDiffusionForceField<DataTypes>::weightVectorPointer::updateResult(size_t degree,typename HighOrderTetrahedralDiffusionForceField<DataTypes>::weightArrayPointer &wap,
																						   typename HighOrderTetrahedralDiffusionForceField<DataTypes>::Vec6 &input) {
	if (degree==2) {
		(*resultQuadratic)+=*(wap.weightArrayQuadratic)*input;
	} else 	if (degree==3) {
		(*resultCubic)+=*(wap.weightArrayCubic)*input;
	} else 	if (degree==4) {
		(*resultQuartic)+=*(wap.weightArrayQuartic)*input;
	} else 	if (degree==5) {
		(*resultQuintic)+=*(wap.weightArrayQuintic)*input;
	}
}
template< class DataTypes>
void HighOrderTetrahedralDiffusionForceField<DataTypes>::weightVectorPointer::updateArray(size_t degree,sofa::helper::vector<Real> &array) {
	Real *ptr;
	if (degree==2) {
		ptr=&(*resultQuadratic)[0];
	} else 	if (degree==3) {
		ptr=&(*resultCubic)[0];
	} else 	if (degree==4) {
		ptr=&(*resultQuartic)[0];
	} else 	if (degree==5) {
		ptr=&(*resultQuintic)[0];
	}
	typename sofa::helper::vector<typename DataTypes::Real>::iterator it=array.begin();
	for(;it!=array.end();it++,ptr++) 
		(*it)+=(*ptr);
}

template< class DataTypes>
void HighOrderTetrahedralDiffusionForceField<DataTypes>::FTCFTetrahedronHandler::applyCreateFunction(unsigned int tetrahedronIndex,
        TetrahedronRestInformation &my_tinfo,
        const Tetrahedron &,
        const sofa::helper::vector<unsigned int> &,
        const sofa::helper::vector<double> &)
{
	if (ff)
	{
		const std::vector< Tetrahedron > &tetrahedronArray=ff->_topology->getTetrahedra() ;
		HighOrderTetrahedronSetTopologyContainer *container=ff->highOrderTetraGeo->getTopologyContainer();
		HighOrderDegreeType degree=container->getDegree();
		size_t nbControlPoints=(degree+1)*(degree+2)*(degree+3)/6;
		size_t nbStiffnessEntries=nbControlPoints*(nbControlPoints-1)/2;
		if (my_tinfo.stiffnessVector.size()!=nbStiffnessEntries) {
			my_tinfo.stiffnessVector.resize(nbStiffnessEntries);
		}
		// set array to zero
		std::fill(my_tinfo.stiffnessVector.begin(),my_tinfo.stiffnessVector.end(),Real());

		//		const std::vector< Edge> &edgeArray=ff->_topology->getEdges() ;
		size_t i,j,k,l,m,n;

		Vec3 point[4];


		Vec6 edgeStiffness;  // the off-diagonal 3x3 block matrices that makes the 12x12 linear elastic matrix

		 const typename HighOrderTetrahedralDiffusionForceField<DataTypes>::MechanicalTypes::VecCoord  &restPosition=ff->mechanicalObject->read(core::ConstVecCoordId::restPosition())->getValue();
		// now computed the stiffness for the HighOrder Tetrahedron
		sofa::helper::vector<TetrahedronIndexVector> tbiArray;

		tbiArray=ff->highOrderTetraGeo->getTopologyContainer()->getTetrahedronIndexArray();

		size_t rank=0;

		///describe the indices of the 4 tetrahedron vertices
		const Tetrahedron &t= tetrahedronArray[tetrahedronIndex];
		//    BaseMeshTopology::EdgesInTetrahedron te=ff->_topology->getEdgesInTetrahedron(tetrahedronIndex);


		// store the point position
		for(j=0; j<4; ++j)
			point[j]=(restPosition)[t[j]];


		if ((ff->integrationMethod==HighOrderTetrahedralDiffusionForceField<DataTypes>::AFFINE_ELEMENT_INTEGRATION) || 
			((ff->d_forceAffineAssemblyForAffineElements.getValue()) && ((ff->highOrderTetraGeo->isBezierTetrahedronAffine(tetrahedronIndex,restPosition  ))))) {


				helper::system::thread::ctime_t startUpdateMat=helper::system::thread::CTime::getTime();
				ff->computeTetrahedronStiffnessEdgeMatrix(point,edgeStiffness);
				if (degree==1) {
				for (rank=0;rank<nbStiffnessEntries;rank++) {
					my_tinfo.stiffnessVector[rank]+=  edgeStiffness[rank];
				}
			} else if (degree==2) {
				Vec45 res=(*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayQuadratic))*edgeStiffness;

				for (rank=0;rank<nbStiffnessEntries;rank++) {
					my_tinfo.stiffnessVector[rank]+=  res[rank];
				}
			} else if (degree==3) {
				Vec190 res=(*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayCubic))*edgeStiffness;

				for (rank=0;rank<nbStiffnessEntries;rank++) {
					my_tinfo.stiffnessVector[rank]+=  res[rank];
				}
			} else if (degree==4) {
				Vec595 res=(*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayQuartic))*edgeStiffness;

				for (rank=0;rank<nbStiffnessEntries;rank++) {
					my_tinfo.stiffnessVector[rank]+=  res[rank];
				}
			} else if (degree==5) {
				Vec1540 res=(*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayQuintic))*edgeStiffness;

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
		} else if (ff->integrationMethod==HighOrderTetrahedralDiffusionForceField<DataTypes>::STANDARD_INTEGRATION){
			sofa::defaulttype::Vec<4,Real> bc;
			helper::system::thread::ctime_t startUpdateMat=helper::system::thread::CTime::getTime();
			// loop through the integration points
			for (i=0;i<ff->numericalIntegrationStiffnessDataArray.size();++i) {

				// the barycentric coordinate
				bc=ff->numericalIntegrationStiffnessDataArray[i].integrationPoint;
				// Compute the local tetrahedron by storing in point the 4 derivatives of the shape functions  
				ff->highOrderTetraGeo->computeNodalValueDerivatives(tetrahedronIndex,bc, restPosition,point);

				Mat3x3 Jacobian,inverseJacobian;
				for (j=0;j<3;++j) {
					for (k=0;k<3;++k) {
						Jacobian[j][k]=point[j][k]-point[3][k];
					}
				}
				invertMatrix(inverseJacobian,Jacobian);
				Real jac=fabs(determinant(Jacobian))*ff->numericalIntegrationStiffnessDataArray[i].integrationWeight;

				helper::vector<Vec3> SVArray;
				for (j=0;j<nbControlPoints;j++) {
					Vec3 sv=inverseJacobian*ff->numericalIntegrationStiffnessDataArray[i].coefficientArray[j];
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

		} else if (ff->integrationMethod==HighOrderTetrahedralDiffusionForceField<DataTypes>::NUMERICAL_INTEGRATION_2){
			helper::system::thread::ctime_t startUpdateMat=helper::system::thread::CTime::getTime();
			sofa::defaulttype::Vec<4,Real> bc;

			// loop through the integration points
			for (i=0;i<ff->numericalIntegrationStiffnessDataArray.size();++i) {

				// the barycentric coordinate
				bc=ff->numericalIntegrationStiffnessDataArray[i].integrationPoint;
				// Compute the local tetrahedron by storing in point the 4 derivatives of the shape functions  
				ff->highOrderTetraGeo->computeNodalValueDerivatives(tetrahedronIndex,bc, restPosition,point);
				// compute the edge stiffness associated with that local tetrahedron
				ff->computeTetrahedronStiffnessEdgeMatrix(point,edgeStiffness);

				// compute the stiffness matrix for all pairs of control points 
				for (rank=0,j=0;j<nbControlPoints;j++) {

					for (k=j+1;k<nbControlPoints;k++,rank++) {

						const Mat4x4  & coeffMatrix=ff->numericalIntegrationStiffnessDataArray[i].weightArray[rank];
						/// add edge stiffness
						for(l=0;l<6;++l) {
							m=edgesInTetrahedronArray[l][0];
							n=edgesInTetrahedronArray[l][1];
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
        }
        else if (ff->integrationMethod == HighOrderTetrahedralDiffusionForceField<DataTypes>::BEZIER_NUMERICAL_INTEGRATION) {
            helper::system::thread::ctime_t startUpdateMat = helper::system::thread::CTime::getTime();
            sofa::defaulttype::Vec<4, Real> bc;

            std::vector< sofa::defaulttype::Vec<6, Real> > reducedStiffness;
            assert(ff->numericalIntegrationStiffnessDataArray.size() > 0);
            size_t numberReducedEntries = ff->numericalIntegrationStiffnessDataArray[0].weightBezierArray.size();
            reducedStiffness.resize(numberReducedEntries);
            // loop through the integration points
            for (i = 0; i < ff->numericalIntegrationStiffnessDataArray.size(); ++i) {

                // the barycentric coordinate
                bc = ff->numericalIntegrationStiffnessDataArray[i].integrationPoint;
                // Compute the local tetrahedron by storing in point the 4 derivatives of the shape functions  
                ff->highOrderTetraGeo->computeNodalValueDerivatives(tetrahedronIndex, bc, restPosition, point);
                // compute the edge stiffness associated with that local tetrahedron
                ff->computeTetrahedronStiffnessEdgeMatrix(point, edgeStiffness);
                for (j = 0; j < numberReducedEntries; ++j) {
                    reducedStiffness[j] += ff->numericalIntegrationStiffnessDataArray[i].weightBezierArray[j] * edgeStiffness;
                }
            }
            size_t r;
            for (i = 0; i < nbStiffnessEntries; ++i) {
                sofa::defaulttype::Vec<16, int>  &mapping = ff->bezierMappingArray[i];
                for (rank = 0, j = 0; j < 4; ++j) {
                    for (k = j+1; k < 4;  ++k,++rank) {
                        r = j * 4 + k;
                        if (mapping[r] >= 0) 
                        my_tinfo.stiffnessVector[i] += ff->bezierCoefficientArray[i][r] * reducedStiffness[(size_t)mapping[r]][rank];
                        r = k * 4 + j;
                        if (mapping[r] >= 0)
                        my_tinfo.stiffnessVector[i] += ff->bezierCoefficientArray[i][r] * reducedStiffness[(size_t)mapping[r]][rank];
                        r = j * 4 + j;
                        if (mapping[r] >= 0)
                        my_tinfo.stiffnessVector[i] -= ff->bezierCoefficientArray[i][r] * reducedStiffness[(size_t)mapping[r]][rank];
                        r = k * 4 + k;
                        if (mapping[r] >= 0)
                        my_tinfo.stiffnessVector[i] -= ff->bezierCoefficientArray[i][r] * reducedStiffness[(size_t)mapping[r]][rank];
                       
                    }
                }
            }
            if (ff->f_printLog.getValue()) {
                helper::system::thread::ctime_t endUpdateMat = helper::system::thread::CTime::getTime();
                ff->totalUpdateMat += endUpdateMat - startUpdateMat;
            }
		} else if (ff->integrationMethod==HighOrderTetrahedralDiffusionForceField<DataTypes>::NUMERICAL_INTEGRATION){
			helper::system::thread::ctime_t startUpdateMat=helper::system::thread::CTime::getTime();
			sofa::defaulttype::Vec<4,Real> bc;

			HighOrderTetrahedralDiffusionForceField<DataTypes>::weightVectorPointer wvp;
			wvp.allocate(degree);
			// loop through the integration points
			for (i=0;i<ff->numericalIntegrationStiffnessDataArray.size();++i) {

				// the barycentric coordinate
				bc=ff->numericalIntegrationStiffnessDataArray[i].integrationPoint;
				// Compute the local tetrahedron by storing in point the 4 derivatives of the shape functions  
				ff->highOrderTetraGeo->computeNodalValueDerivatives(tetrahedronIndex,bc, restPosition,point);
				// compute the edge stiffness associated with that local tetrahedron
				ff->computeTetrahedronStiffnessEdgeMatrix(point,edgeStiffness);
				wvp.updateResult(degree,ff->numericalIntegrationStiffnessDataArray[i].arrayPointer,edgeStiffness);

			}
			wvp.updateArray(degree,my_tinfo.stiffnessVector);
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
					std::cerr<< "Stiffness entry ["<<(unsigned int)tbiArray[l][0]<<" "<< (unsigned int)tbiArray[l][1]<<" "<< (unsigned int)tbiArray[l][2]<<" "<<(unsigned int) tbiArray[l][3]<< "]["<<
						(unsigned int)tbiArray[m][0]<<" "<< (unsigned int)tbiArray[m][1]<<" "<< (unsigned int)tbiArray[m][2]<<" "<< (unsigned int)tbiArray[m][3]<<"]="<<my_tinfo.stiffnessVector[rank]<<std::endl;

				}
			}
		}
#endif
	}
}



template <class DataTypes> HighOrderTetrahedralDiffusionForceField<DataTypes>::HighOrderTetrahedralDiffusionForceField()
    : tetrahedronInfo(initData(&tetrahedronInfo, "tetrahedronInfo", "Internal tetrahedron data"))
    , updateMatrix(true)
    , d_diffusivity(initData(&d_diffusivity,(Real)1000.,"diffusivity","diffusivity for isotropic diffusion"))
	, d_anisotropy(initData(&d_anisotropy,std::string("isotropy"),"anisotropy","\"isotropy\" or \"transverseIsotropy\" or \"orthotropy\" as an anisotropy of diffusion"))
	, numericalIntegrationOrder( initData(&numericalIntegrationOrder,(size_t)2,"integrationOrder","The order of integration for numerical integration"))
	, d_integrationMethod( initData(&d_integrationMethod,std::string("analytical"),"integrationMethod","\"analytical\" if closed form expression for affine element, \"numerical\" if numerical integration is chosen,  \"standard\" if standard integration is chosen"))
	, d_anisotropyParameter( initData(&d_anisotropyParameter,ParameterArray(),"anisotropyParameters","diffusivity parameters if the diffusion is anisotropic "))
	, d_anisotropyDirection( initData(&d_anisotropyDirection,AnisotropyDirectionArray(),"anisotropyDirection","Directions of anisotropy"))
	, numericalIntegrationMethod( initData(&numericalIntegrationMethod, std::string("Tetrahedron Gauss"),"numericalIntegrationMethod","The type of numerical integration method chosen"))
    , d_assemblyTime(initData(&d_assemblyTime,(Real)0,"assemblyTime","the time spent in assembling the stiffness matrix. Only updated if printLog is set to true"))
    , d_forceAffineAssemblyForAffineElements(initData(&d_forceAffineAssemblyForAffineElements,true,"forceAffineAssemblyForAffineElements","if true affine tetrahedra are always assembled with the closed form formula, Otherwise use the method defined in integrationMethod"))
	, tetrahedronHandler(NULL)
{
    tetrahedronHandler = new FTCFTetrahedronHandler(this,&tetrahedronInfo);
}

template <class DataTypes> HighOrderTetrahedralDiffusionForceField<DataTypes>::~HighOrderTetrahedralDiffusionForceField()
{
    if (tetrahedronHandler) delete tetrahedronHandler;
}


template <class DataTypes> void HighOrderTetrahedralDiffusionForceField<DataTypes>::init()
{
    //	serr << "initializing HighOrderTetrahedralDiffusionForceField" << sendl;
    this->Inherited::init();

    _topology = this->getContext()->getMeshTopology();
	this->getContext()->get(highOrderTetraGeo);

    if ((_topology->getNbTetrahedra()==0) || (!highOrderTetraGeo))
    {
        serr << "ERROR(HighOrderTetrahedralDiffusionForceField): object must have a Tetrahedral Set Topology and a HighOrderTetrahedronSetGeometryAlgorithms component "<<sendl;
        return;
    }
   
	if (d_integrationMethod.getValue() == "analytical")
        integrationMethod= AFFINE_ELEMENT_INTEGRATION;
    else if (d_integrationMethod.getValue() == "numerical") 
		integrationMethod= NUMERICAL_INTEGRATION;
	else if (d_integrationMethod.getValue() == "numerical2") 
		integrationMethod= NUMERICAL_INTEGRATION_2;
    else if (d_integrationMethod.getValue() == "bezierNumerical")
        integrationMethod = BEZIER_NUMERICAL_INTEGRATION;
	else if (d_integrationMethod.getValue() == "standard") 
		integrationMethod= STANDARD_INTEGRATION;
    else
    {
        serr << "cannot recognize method "<< d_integrationMethod.getValue() << ". Must be either \"analytical\" or \"numerical\" or \"bezierNumerical\" or \"standard\"" << sendl;
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

    std::set<typename topology::NumericalIntegrationDescriptor<Real, 3>::QuadratureMethod> qmSet = highOrderTetraGeo->getTetrahedronNumericalIntegrationDescriptor().getQuadratureMethods();
    if (qmSet.count(numericalIntegrationMethod.getValue()) == 0) {
        serr << "cannot recognize numerical integration method  " << numericalIntegrationMethod.getValue() << sendl;
    }

    helper::vector<TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());
    tetrahedronInf.resize(_topology->getNbTetrahedra());


	/// Get the mechanical object containing the mesh position in 3D
	sofa::core::objectmodel::Tag mechanicalTag(m_tagMeshMechanics.getValue());
	this->getContext()->get(mechanicalObject, mechanicalTag,sofa::core::objectmodel::BaseContext::SearchUp);
	if (mechanicalObject==NULL)
	{
		serr<<"ERROR(HighOrderTetrahedralDiffusionForceField): cannot find the mechanical object."<<sendl;
		std::cout<<"mechanicalObj = "<<mechanicalObject<<std::endl;
		return;
	}

	helper::system::thread::ctime_t startAssembly=helper::system::thread::CTime::getTime();
	totalUpdateMat=0;

	size_t i;

	topology::HighOrderDegreeType degree=highOrderTetraGeo->getTopologyContainer()->getDegree();
	size_t nbControlPoints=(degree+1)*(degree+2)*(degree+3)/6;
	sofa::helper::vector<topology::TetrahedronIndexVector> tbiArray;
	tbiArray=highOrderTetraGeo->getTopologyContainer()->getTetrahedronIndexArray();
	topology::TetrahedronIndexVector tbi1,tbi2;
	
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
			std::vector<Real> coeffArray(6);
			Mat4x4 coeffMatrix;
			size_t j,k,l,m,n,rank;
			for (rank=0,j=0;j<nbControlPoints;j++) {
				tbi1=tbiArray[j];
				for (k=j+1;k<nbControlPoints;k++,rank++) {
					tbi2=tbiArray[k];
					coeffMatrix=highOrderTetraGeo->getAffineStiffnessCoefficientMatrix(tbi1,tbi2);
					// substract the diagonal terms such that only edge stiffness are used
					for(l=0; l<4; ++l){
						for(m=l+1; m<4; ++m){
							coeffMatrix[l][m]+= coeffMatrix[m][l]-(coeffMatrix[l][l]+coeffMatrix[m][m]);						

						}
					}
					if (degree==2) {
						for(l=0; l<6; ++l){
							m=edgesInTetrahedronArray[l][0];
							n=edgesInTetrahedronArray[l][1];

							(*(affineStiffnessCoefficientPreStoredArray.weightArrayQuadratic))[rank][l]=coeffMatrix[m][n];

						}
					} else if (degree==3) {
						for(l=0; l<6; ++l){
							m=edgesInTetrahedronArray[l][0];
							n=edgesInTetrahedronArray[l][1];

							(*(affineStiffnessCoefficientPreStoredArray.weightArrayCubic))[rank][l]=coeffMatrix[m][n];

						}
					} else if (degree==4) {
						for(l=0; l<6; ++l){
							m=edgesInTetrahedronArray[l][0];
							n=edgesInTetrahedronArray[l][1];

							(*(affineStiffnessCoefficientPreStoredArray.weightArrayQuartic))[rank][l]=coeffMatrix[m][n];

						}
					} else if (degree==5) {
						for(l=0; l<6; ++l){
							m=edgesInTetrahedronArray[l][0];
							n=edgesInTetrahedronArray[l][1];

							(*(affineStiffnessCoefficientPreStoredArray.weightArrayQuintic))[rank][l]=coeffMatrix[m][n];

						}

					} else {
						Vec6 coeff;
						for(l=0; l<6; ++l){
							m=edgesInTetrahedronArray[l][0];
							n=edgesInTetrahedronArray[l][1];
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
			topology::NumericalIntegrationDescriptor<Real,4> &nid=highOrderTetraGeo->getTetrahedronNumericalIntegrationDescriptor();
			typename topology::NumericalIntegrationDescriptor<Real,4>::QuadraturePointArray qpa=nid.getQuadratureMethod((typename topology::NumericalIntegrationDescriptor<Real,4>::QuadratureMethod)numericalIntegrationMethod.getValue(),
				numericalIntegrationOrder.getValue());
			size_t i,j,k,l,m;
			sofa::defaulttype::Vec<4,Real> bc;
			Real weight;
			Mat4x4 coeffMatrix;

			// loop through the integration points
			for (i=0;i<qpa.size();++i) {
				NumericalIntegrationStiffnessData nimd;
				typename topology::NumericalIntegrationDescriptor<Real,4>::QuadraturePoint qp=qpa[i];
				// the barycentric coordinate
				nimd.integrationPoint=qp.first;
				// the weight of the integration point
				weight=qp.second;
				nimd.integrationWeight=qp.second;
				if ((integrationMethod== NUMERICAL_INTEGRATION) && (degree<6)) {
					nimd.arrayPointer.allocate(degree);
				}
				size_t n;

				std::vector<Vec4> shapeFunctionDerivativeArray;
				for(j=0;j<tbiArray.size();++j) {
					Vec4 deriv=highOrderTetraGeo->computeShapeFunctionDerivatives(tbiArray[j],qp.first);
					shapeFunctionDerivativeArray.push_back(deriv);
					Vec3 der(deriv[0]-deriv[3],deriv[1]-deriv[3],deriv[2]-deriv[3]);
					nimd.coefficientArray.push_back(der);
				}
				size_t rank;
				for(rank=0,j=0;j<tbiArray.size();++j) {
					for(k=j+1;k<tbiArray.size();++k,++rank) {
						coeffMatrix=dyad(shapeFunctionDerivativeArray[j],shapeFunctionDerivativeArray[k])*6*weight;
						for(l=0; l<4; ++l){
							for(m=l+1; m<4; ++m){
								coeffMatrix[l][m]+= coeffMatrix[m][l]-(coeffMatrix[l][l]+coeffMatrix[m][m]);
							}
						}
						if (integrationMethod== NUMERICAL_INTEGRATION) { 
							if (degree>5)  {
								Vec6 coeffVec;
								for(l=0; l<6; ++l){
									m=edgesInTetrahedronArray[l][0];
									n=edgesInTetrahedronArray[l][1];
									coeffVec[l]=coeffMatrix[m][n];								
								}
								nimd.weightVectorizedArray.push_back(coeffVec);

							} else {
								if (degree==2) {
									//							nimd.weightArrayQuadratic=new Mat45x6;
									for(l=0; l<6; ++l){
										m=edgesInTetrahedronArray[l][0];
										n=edgesInTetrahedronArray[l][1];
										(*(nimd.arrayPointer.weightArrayQuadratic))[rank][l]=coeffMatrix[m][n];
									}
								} else 	if (degree==3) {
									//							nimd.weightArrayCubic=new Mat45x6;
									for(l=0; l<6; ++l){
										m=edgesInTetrahedronArray[l][0];
										n=edgesInTetrahedronArray[l][1];
										(*(nimd.arrayPointer.weightArrayCubic))[rank][l]=coeffMatrix[m][n];
									}
								} else 	if (degree==4) {

									//							nimd.weightArrayQuartic=new Mat45x6;
									for(l=0; l<6; ++l){
										m=edgesInTetrahedronArray[l][0];
										n=edgesInTetrahedronArray[l][1];
										(*(nimd.arrayPointer.weightArrayQuartic))[rank][l]=coeffMatrix[m][n];
									}
								} else 	if (degree==5) {
									//							nimd.weightArrayQuintic=new Mat45x6;
									for(l=0; l<6; ++l){
										m=edgesInTetrahedronArray[l][0];
										n=edgesInTetrahedronArray[l][1];
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
        if (integrationMethod == BEZIER_NUMERICAL_INTEGRATION)
        {
            /// first fill the first vector independent from the integration points
            BezierTetrahedronSetGeometryAlgorithms<MechanicalTypes> *bezierTetraGeo = dynamic_cast<BezierTetrahedronSetGeometryAlgorithms<MechanicalTypes> *> (highOrderTetraGeo);
            if (bezierTetraGeo == NULL) {
                serr << "Could not find any BezierTetrahedronSetGeometryAlgorithms while using BEZIER_NUMERICAL_INTEGRATION" << sendl;
                return;
            }
            bezierCoefficientArray.clear();
            /// store the coefficient that are independent from the integration point
            topology::TetrahedronIndexVector tbi1Copy, tbi2Copy;
            size_t rank,r;
            size_t i, j, k, l, m;
            for (rank = 0, j = 0; j < tbiArray.size(); ++j) {
                for (k = j + 1; k < tbiArray.size(); ++k, ++rank) {
                    Vec16 weight;
                    tbi1 = tbiArray[j];
                    tbi2=  tbiArray[k];
                    for (r=0,l = 0; l < 4; ++l) {
                        for (m = 0; m < 4; ++m,++r) {
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
            sofa::helper::vector<topology::TetrahedronIndexVector> tbiArray2;
            tbiArray2 = bezierTetraGeo->getTopologyContainer()->getTetrahedronIndexArrayOfGivenDegree(2 * degree - 2);

            // first create a map to speed up the assignment of index from  TetrahedronIndexVector
            std::map<topology::TetrahedronIndexVector, size_t> tivMap;
            std::map<topology::TetrahedronIndexVector, size_t>::iterator itmap;
            for ( j = 0; j < tbiArray2.size(); ++j) {
                tivMap.insert(std::make_pair(tbiArray2[j],j));
            }
            /// now fills the bezierMappingArray
            for (rank = 0, j = 0; j < tbiArray.size(); ++j) {
                for (k = j + 1; k < tbiArray.size(); ++k, ++rank) {
                    tbi1= tbiArray[j] + tbiArray[k];
                    Vec16Int mapping;
                    for (r=0,l = 0; l<4; ++l) {
                        for (m = 0; m<4; ++m,++r) {
                            
                            if (((tbi1[l]*tbi1[m]) == 0) || (((l == m) && (tbi1[l] < 2)))) {
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
            topology::NumericalIntegrationDescriptor<Real, 4> &nid = highOrderTetraGeo->getTetrahedronNumericalIntegrationDescriptor();
            typename topology::NumericalIntegrationDescriptor<Real, 4>::QuadraturePointArray qpa = nid.getQuadratureMethod((typename topology::NumericalIntegrationDescriptor<Real, 4>::QuadratureMethod)numericalIntegrationMethod.getValue(),
                numericalIntegrationOrder.getValue());
          
            sofa::defaulttype::Vec<4, Real> bc;
            Real weight,fac;
           
            fac = (Real)lfactorial(degree - 1)*(Real)lfactorial(degree - 1) / (Real)lfactorial(2 * degree - 2);

            // loop through the integration points
            for (i = 0; i<qpa.size(); ++i) {
                NumericalIntegrationStiffnessData nimd;
                typename topology::NumericalIntegrationDescriptor<Real, 4>::QuadraturePoint qp = qpa[i];
                // the barycentric coordinate
                nimd.integrationPoint = qp.first;
                // the weight of the integration point
                weight = qp.second;
                nimd.integrationWeight = qp.second;
  
                nimd.weightBezierArray.resize(tbiArray2.size());
                for (j = 0; j < tbiArray2.size(); ++j) {
                    nimd.weightBezierArray[j] = 6*degree*degree*fac*bezierTetraGeo->computeShapeFunctionOfGivenDegree(tbiArray2[j], qp.first, 2 * degree - 2);
                    nimd.weightBezierArray[j] *= weight*factorialTVI(tbiArray2[j]);
                }

                numericalIntegrationStiffnessDataArray.push_back(nimd);
            }
        }
		if (integrationMethod== STANDARD_INTEGRATION) 
		{
			numericalIntegrationStiffnessDataArray.clear();
			/// get value of integration points0
			topology::NumericalIntegrationDescriptor<Real,4> &nid=highOrderTetraGeo->getTetrahedronNumericalIntegrationDescriptor();
			typename topology::NumericalIntegrationDescriptor<Real,4>::QuadraturePointArray qpa=nid.getQuadratureMethod((typename topology::NumericalIntegrationDescriptor<Real,4>::QuadratureMethod)numericalIntegrationMethod.getValue(),
				numericalIntegrationOrder.getValue());
			size_t i,j,k;
			sofa::defaulttype::Vec<4,Real> bc;
			Real weight;
			Mat4x4 coeffMatrix;

			// loop through the integration points
			for (i=0;i<qpa.size();++i) {
				NumericalIntegrationStiffnessData nimd;
				typename topology::NumericalIntegrationDescriptor<Real,4>::QuadraturePoint qp=qpa[i];
				// the barycentric coordinate
				nimd.integrationPoint=qp.first;
				// the weight of the integration point	
				nimd.integrationWeight=qp.second;

				std::vector<Vec4> shapeFunctionDerivativeArray;
				for(j=0;j<tbiArray.size();++j) {
					Vec4 deriv=highOrderTetraGeo->computeShapeFunctionDerivatives(tbiArray[j],qp.first);
					Vec3 der(deriv[0]-deriv[3],deriv[1]-deriv[3],deriv[2]-deriv[3]);
					nimd.coefficientArray.push_back(der);
				}
				numericalIntegrationStiffnessDataArray.push_back(nimd);
			}



		}
	}
	    /// initialize the data structure associated with each tetrahedron
    for (i=0; i<_topology->getNbTetrahedra(); ++i)
    {
        tetrahedronHandler->applyCreateFunction(i,tetrahedronInf[i],_topology->getTetrahedron(i),
                (const helper::vector< unsigned int > )0,
                (const helper::vector< double >)0);
    }
	if (this->f_printLog.getValue()) {
		helper::system::thread::ctime_t endAssembly=helper::system::thread::CTime::getTime();
		std::cerr<< "Assembly time="<< ((endAssembly-startAssembly)/(Real)helper::system::thread::CTime::getRefTicksPerSec()) << std::endl;
		std::cerr<<" total update mat="<<((totalUpdateMat)/(Real)helper::system::thread::CTime::getRefTicksPerSec()) << std::endl;
		d_assemblyTime.setValue(((endAssembly-startAssembly)/(Real)helper::system::thread::CTime::getRefTicksPerSec()));
	}
    /// set the call back function upon creation of a tetrahedron
    tetrahedronInfo.createTopologicalEngine(_topology,tetrahedronHandler);
    tetrahedronInfo.registerTopologicalData();
    tetrahedronInfo.endEdit();

	updateTopologyInfo=true;

}


template <class DataTypes>
void HighOrderTetrahedralDiffusionForceField<DataTypes>::updateTopologyInformation()
{
    int i;
    unsigned int j;

    int nbTetrahedra=_topology->getNbTetrahedra();

    TetrahedronRestInformation *tetinfo;

    helper::vector<typename HighOrderTetrahedralDiffusionForceField<DataTypes>::TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());
   


    for(i=0; i<nbTetrahedra; i++ )
    {
        tetinfo=&tetrahedronInf[i];
        /// describe the jth edge index of triangle no i
        const EdgesInTetrahedron &tea= _topology->getEdgesInTetrahedron(i);
        /// describe the jth vertex index of triangle no i
        const Tetrahedron &ta= _topology->getTetrahedron(i);

        for (j=0; j<4; ++j)
        {
            tetinfo->v[j]=ta[j];
        }


    }
    updateTopologyInfo=false;
    tetrahedronInfo.endEdit();
}
template<class DataTypes>
void HighOrderTetrahedralDiffusionForceField<DataTypes>::computeDiffusivityTensor() 											
{
	const helper::vector<Real> & anisotropyParameter=d_anisotropyParameter.getValue();
	const helper::vector<Vec3> & anisotropyDirection=d_anisotropyDirection.getValue();
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
void HighOrderTetrahedralDiffusionForceField<DataTypes>::computeTetrahedronStiffnessEdgeMatrix(const Vec3 point[4],
																							   Vec6 &edgeStiffness)
{
	Vec3 shapeVector[4];
	/// compute 6 times the rest volume
	Real volume=dot(cross(point[1]-point[0],point[2]-point[0]),point[0]-point[3]);
	/// store the rest volume
//	my_tinfo.restVolume=volume/6;

	size_t j,k,l;
	// store shape vectors at the rest configuration
	for(j=0; j<4; ++j)
	{
		if ((j%2)==0)
			shapeVector[j]=cross(point[(j+2)%4] - point[(j+1)%4],point[(j+3)%4] - point[(j+1)%4])/volume;
		else
			shapeVector[j]= -cross(point[(j+2)%4] - point[(j+1)%4],point[(j+3)%4] - point[(j+1)%4])/volume;
		
	}
	if (diffusionSymmetry==ISOTROPIC) {
		Real diff=d_diffusivity.getValue()*fabs(volume)/6;


		/// compute the edge stiffness of the linear elastic material
		for(j=0; j<6; ++j)
		{
			k=edgesInTetrahedronArray[j][0];
			l=edgesInTetrahedronArray[j][1];
			// the linear stiffness matrix using shape vectors and Lame coefficients
			edgeStiffness[j]=diff*dot(shapeVector[l],shapeVector[k]);
		}
	} else {
		Mat3x3 diff=diffusionTensor*fabs(volume)/6;
		for(j=0; j<6; ++j)
		{
			k=edgesInTetrahedronArray[j][0];
			l=edgesInTetrahedronArray[j][1];
			edgeStiffness[j]=dot(shapeVector[l],diff*shapeVector[k]);
		}
	}
}

template <class DataTypes>
const helper::vector<typename HighOrderTetrahedralDiffusionForceField<DataTypes>::Real> & 
	HighOrderTetrahedralDiffusionForceField<DataTypes>::getStiffnessArray(
	const size_t i,
	const typename HighOrderTetrahedralDiffusionForceField<DataTypes>::TetrahedronRestInformation *restTetra)
{
	return(restTetra->stiffnessVector);
}

template <class DataTypes>
void HighOrderTetrahedralDiffusionForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* mparams /* PARAMS FIRST */, 
																	 DataVecDeriv &  dataF, const DataVecCoord &  dataX , const DataVecDeriv & /*dataV*/ )
{

    VecDeriv& f        = *(dataF.beginEdit());
    const VecCoord& x  =   dataX.getValue()  ;
	
    size_t i,j,k,l,v0,v1,rank,p;
	size_t nbTetrahedra=_topology->getNbTetrahedra();
	HighOrderDegreeType degree=highOrderTetraGeo->getTopologyContainer()->getDegree();
	size_t nbControlPoints=(degree+1)*(degree+2)*(degree+3)/6;

	HighOrderTetrahedronSetTopologyContainer::VecPointID indexArray;


    if (updateTopologyInfo)
    {
        updateTopologyInformation();
    }
    helper::vector<TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());
    TetrahedronRestInformation *tetinfo;

	
	
	Coord dpos,sv;
		
	for(i=0; i<nbTetrahedra; i++ )
	{
		tetinfo=&tetrahedronInf[i];
		const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=
			highOrderTetraGeo->getTopologyContainer()->getGlobalIndexArrayOfControlPoints(i);
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
    tetrahedronInfo.endEdit();

    dataF.endEdit();

}


template <class DataTypes>
void HighOrderTetrahedralDiffusionForceField<DataTypes>::addDForce(const sofa::core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv&   datadF , const DataVecDeriv&   datadX )
{
    VecDeriv& df       = *(datadF.beginEdit());
    const VecCoord& dx =   datadX.getValue()  ;
    Real kFactor = (Real)mparams->kFactor();
    const helper::vector<TetrahedronRestInformation>& tetrahedronInf = tetrahedronInfo.getValue();
	Coord dpos;
    size_t i,j,k,v0,v1,rank;
    size_t nbTetrahedra=_topology->getNbTetrahedra();
    
	HighOrderDegreeType degree=highOrderTetraGeo->getTopologyContainer()->getDegree();
	size_t nbControlPoints=(degree+1)*(degree+2)*(degree+3)/6;
	
    const TetrahedronRestInformation *tetinfo;
	

	for(i=0; i<nbTetrahedra; i++ )
	{
		tetinfo=&tetrahedronInf[i];
		const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=
			highOrderTetraGeo->getTopologyContainer()->getGlobalIndexArrayOfControlPoints(i);
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
SReal HighOrderTetrahedralDiffusionForceField<DataTypes>::getPotentialEnergy(const core::MechanicalParams*, const DataVecCoord&) const
{
    serr << "ERROR("<<this->getClassName()<<"): getPotentialEnergy( const MechanicalParams*, const DataVecCoord& ) not implemented." << sendl;
    return 0.0;
}




template<class DataTypes>
void HighOrderTetrahedralDiffusionForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
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

#endif // SOFA_COMPONENT_FORCEFIELD_HighOrderTetrahedralDiffusionForceField_INL
