
#ifndef SOFA_COMPONENT_FORCEFIELD_BEZIERTETRAHEDRALCOROTATIONALFEMFORCEFIELD_INL
#define SOFA_COMPONENT_FORCEFIELD_BEZIERTETRAHEDRALCOROTATIONALFEMFORCEFIELD_INL

#include "BezierTetrahedralCorotationalFEMForceField.h"
#include <sofa/core/visual/VisualParams.h>
#include <fstream> // for reading the file
#include <iostream> //for debugging
#include <sofa/helper/gl/template.h>
#include <SofaBaseTopology/TopologyData.inl>
#include <HighOrderTetrahedronSetGeometryAlgorithms.h>
#include <sofa/core/behavior/ForceField.inl>
#include <SofaBaseTopology/CommonAlgorithms.h>
#include <sofa/helper/decompose.h>

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
const unsigned int myedgesInTetrahedronArray[6][2] = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}};



template< class DataTypes>
void BezierTetrahedralCorotationalFEMForceField<DataTypes>::FTCFTetrahedronHandler::applyCreateFunction(unsigned int tetrahedronIndex,
        TetrahedronRestInformation &my_tinfo,
        const Tetrahedron &,
        const sofa::helper::vector<unsigned int> &,
        const sofa::helper::vector<double> &)
{
	if (ff)
	{
		const std::vector< Tetrahedron > &tetrahedronArray=ff->_topology->getTetrahedra() ;
		HighOrderTetrahedronSetTopologyContainer *container=ff->bezierTetraGeo->getTopologyContainer();
		HighOrderDegreeType degree=container->getDegree();
		size_t nbControlPoints=(degree+1)*(degree+2)*(degree+3)/6;
		size_t nbStiffnessEntries=nbControlPoints*(nbControlPoints-1)/2;
		if (my_tinfo.stiffnessVector.size()!=nbStiffnessEntries) {
			my_tinfo.stiffnessVector.resize(nbStiffnessEntries);
		}
		// set array to zero
		std::fill(my_tinfo.stiffnessVector.begin(),my_tinfo.stiffnessVector.end(),Mat3x3());

		//		const std::vector< Edge> &edgeArray=ff->_topology->getEdges() ;
		size_t i,j,k,l,m,n;

		typename DataTypes::Coord point[4];
		typedef typename Mat<3,3, DataTypes::Real> Mat3x3;
		typedef typename Mat<3,6, DataTypes::Real> Mat3x6;

		Mat3x3 edgeStiffness[6];  // the off-diagonal 3x3 block matrices that makes the 12x12 linear elastic matrix

		const typename DataTypes::VecCoord restPosition=ff->mstate->read(core::ConstVecCoordId::restPosition())->getValue();
		// now computed the stiffness for the Bezier Tetrahedron
		sofa::helper::vector<TetrahedronIndexVector> tbiArray;

		tbiArray=ff->bezierTetraGeo->getTopologyContainer()->getTetrahedronIndexArray();

		size_t rank=0;
		if (ff->d_oneRotationPerIntegrationPoint.getValue()) {
			// one rotation per integration point
			sofa::defaulttype::Vec<4,Real> bc;
			my_tinfo.integrationPointsStiffnessVector.resize(ff->numericalIntegrationStiffnessDataArray.size());
			my_tinfo.integrationPointsRestEdgeVector.resize(ff->numericalIntegrationStiffnessDataArray.size());
			my_tinfo.integrationPointsRestRotationArray.resize(ff->numericalIntegrationStiffnessDataArray.size());

			// loop through the integration points
			for (i=0;i<ff->numericalIntegrationStiffnessDataArray.size();++i) {

				// the barycentric coordinate
				bc=ff->numericalIntegrationStiffnessDataArray[i].integrationPoint;
				// Compute the local tetrahedron by storing in point the 4 derivatives of the shape functions  
				ff->bezierTetraGeo->computeNodalValueDerivatives(tetrahedronIndex,bc, restPosition,point);

				// initialize rotation of element
				if (ff->decompositionMethod==QR_DECOMPOSITION) {
					helper::vector<Coord> restEdgeVector(6);				
					Coord restEdgeVectorTmp[6];

					for(j=0; j<6; ++j){
						k=edgesInTetrahedronArray[j][0];
						l=edgesInTetrahedronArray[j][1];

						// store the rest edge vector
						restEdgeVector[j]=point[l]-point[k];
						restEdgeVectorTmp[j]=restEdgeVector[j];
					}
					my_tinfo.integrationPointsRestEdgeVector[i]=restEdgeVector;

					// compute the rotation matrix of the initial tetrahedron for the QR decomposition
					computeQRRotation(my_tinfo.integrationPointsRestRotationArray[i],
						restEdgeVectorTmp);	
				} else 	if (ff->decompositionMethod==POLAR_DECOMPOSITION_MODIFIED) {
					Mat3x3 Transformation;
					Transformation[0]=point[1]-point[0];
					Transformation[1]=point[2]-point[0];
					Transformation[2]=point[3]-point[0];
					helper::Decompose<Real>::polarDecomposition( Transformation, my_tinfo.integrationPointsRestRotationArray[i] );
				}


				// compute the edge stiffness associated with that local tetrahedron
				ff->computeTetrahedronStiffnessEdgeMatrix(point,edgeStiffness);
				helper::vector<Mat3x3> stiffVector(6);
				for(j=0; j<6; ++j){
					stiffVector[j]=edgeStiffness[j];
				}
				my_tinfo.integrationPointsStiffnessVector[i]=stiffVector;
				my_tinfo.rotatedStiffnessVector.resize(nbControlPoints*(nbControlPoints-1)/2);
			}

		} else {
			///describe the indices of the 4 tetrahedron vertices
			const Tetrahedron &t= tetrahedronArray[tetrahedronIndex];
			//    BaseMeshTopology::EdgesInTetrahedron te=ff->_topology->getEdgesInTetrahedron(tetrahedronIndex);


			// store the point position
			for(j=0; j<4; ++j)
				point[j]=(restPosition)[t[j]];

			if (ff->decompositionMethod==QR_DECOMPOSITION) {
				for(j=0; j<6; ++j){
					k=edgesInTetrahedronArray[j][0];
					l=edgesInTetrahedronArray[j][1];

					// store the rest edge vector
					my_tinfo.restEdgeVector[j]=point[l]-point[k];
				}
				// compute the rotation matrix of the initial tetrahedron for the QR decomposition
				computeQRRotation(my_tinfo.restRotation,my_tinfo.restEdgeVector);
			} else 	if (ff->decompositionMethod==POLAR_DECOMPOSITION_MODIFIED) {
				Mat3x3 Transformation;
				Transformation[0]=point[1]-point[0];
				Transformation[1]=point[2]-point[0];
				Transformation[2]=point[3]-point[0];
				helper::Decompose<Real>::polarDecomposition( Transformation, my_tinfo.restRotation );
			}

			if (ff->integrationMethod==BezierTetrahedralCorotationalFEMForceField<DataTypes>::AFFINE_ELEMENT_INTEGRATION) {

				ff->computeTetrahedronStiffnessEdgeMatrix(point,edgeStiffness);

				for (rank=0,j=0;j<nbControlPoints;j++) {
					for (k=j+1;k<nbControlPoints;k++,rank++) {	
						Mat3x3 stiff;
						const Mat4x4 & wMat=ff->affineStiffnessCoefficientArray[rank];
						/// add edge stiffness
						for(l=0;l<6;++l) {
							m=edgesInTetrahedronArray[l][0];
							n=edgesInTetrahedronArray[l][1];
							if (wMat[m][n]!=0)
								my_tinfo.stiffnessVector[rank]+= edgeStiffness[l]*wMat[m][n];
							if (wMat[n][m]!=0)
								my_tinfo.stiffnessVector[rank]+=edgeStiffness[l].transposed()*wMat[n][m];
						}

					}
				}
				
			} else if (ff->integrationMethod==BezierTetrahedralCorotationalFEMForceField<DataTypes>::STANDARD_INTEGRATION){
				sofa::defaulttype::Vec<4,Real> bc;

				// loop through the integration points
				for (i=0;i<ff->numericalIntegrationStiffnessDataArray.size();++i) {

					// the barycentric coordinate
					bc=ff->numericalIntegrationStiffnessDataArray[i].integrationPoint;
					// Compute the local tetrahedron by storing in point the 4 derivatives of the shape functions  
					ff->bezierTetraGeo->computeNodalValueDerivatives(tetrahedronIndex,bc, restPosition,point);

					Mat3x3 Jacobian,inverseJacobian;
					for (j=0;j<3;++j) {
						for (k=0;k<3;++k) {
							Jacobian[j][k]=point[j][k]-point[3][k];
						}
					}
					invertMatrix(inverseJacobian,Jacobian);
					Real jac=fabs(determinant(Jacobian))*ff->numericalIntegrationStiffnessDataArray[i].integrationWeight;

					helper::vector<Mat6x3> SDArray;
					for (j=0;j<nbControlPoints;j++) {
						Coord sv=inverseJacobian*ff->numericalIntegrationStiffnessDataArray[i].coefficientArray[j];
						Mat6x3 strainDisplacement;
						strainDisplacement[0][0]=sv[0];strainDisplacement[1][1]=sv[1];strainDisplacement[2][2]=sv[2];
						strainDisplacement[3][1]=sv[0];strainDisplacement[3][0]=sv[1];strainDisplacement[4][1]=sv[2];
						strainDisplacement[5][2]=sv[0];strainDisplacement[4][2]=sv[1];strainDisplacement[5][0]=sv[2];
						SDArray.push_back(strainDisplacement);	
					}

					for (rank=0,j=0;j<nbControlPoints;j++) {

						for (k=j+1;k<nbControlPoints;k++,rank++) {
							Mat3x3 edgeStiffness=(SDArray[j].transposed()*((ff->elasticityTensor)*SDArray[k]));
							edgeStiffness*= jac;
							my_tinfo.stiffnessVector[rank]+=edgeStiffness.transposed();
						}
					}
				}

			} else if (ff->integrationMethod==BezierTetrahedralCorotationalFEMForceField<DataTypes>::NUMERICAL_INTEGRATION){

				sofa::defaulttype::Vec<4,Real> bc;



				// loop through the integration points
				for (i=0;i<ff->numericalIntegrationStiffnessDataArray.size();++i) {

					// the barycentric coordinate
					bc=ff->numericalIntegrationStiffnessDataArray[i].integrationPoint;
					// Compute the local tetrahedron by storing in point the 4 derivatives of the shape functions  
					ff->bezierTetraGeo->computeNodalValueDerivatives(tetrahedronIndex,bc, restPosition,point);
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
								if (coeffMatrix[n][m]!=0) {
									my_tinfo.stiffnessVector[rank]+= edgeStiffness[l].transposed()*coeffMatrix[n][m];
								}
									
							}
						}
					}



				}
				//			if (tetrahedronIndex==10) {
				//			std::cerr<<"stiffness at rank 4= "<<my_tinfo.stiffnessVector[4]<<std::endl;
				//		}

			}		
#ifdef _DEBUG
			if (ff->f_printLog.getValue()) {
				std::cerr<< " lambda="<<ff->getLambda() << " mu="<<ff->getMu() << std::endl;


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

}

template <class DataTypes> BezierTetrahedralCorotationalFEMForceField<DataTypes>::BezierTetrahedralCorotationalFEMForceField()
    : tetrahedronInfo(initData(&tetrahedronInfo, "tetrahedronInfo", "Internal tetrahedron data"))
    , _initialPoints(0)
    , updateMatrix(true)
    , d_method(initData(&d_method,std::string("linear"),"method","method for rotation computation :\"qr\" (by QR) or \"polar\" or \"polar2\" or \"none\" (Linear elastic)"))
    , d_poissonRatio(initData(&d_poissonRatio,(Real)0.3,"poissonRatio","Poisson ratio in Hooke's law"))
    , d_youngModulus(initData(&d_youngModulus,(Real)1000.,"youngModulus","Young modulus in Hooke's law"))
	, numericalIntegrationOrder( initData(&numericalIntegrationOrder,(size_t)2,"integrationOrder","The order of integration for numerical integration"))
	, d_integrationMethod( initData(&d_integrationMethod,std::string("analytical"),"integrationMethod","\"analytical\" if closed form expression for affine element, \"numerical\" if numerical integration is chosen,  \"standard\" if standard integration is chosen"))
	, numericalIntegrationMethod( initData(&numericalIntegrationMethod,(size_t)0,"numericalIntegrationMethod","The type of numerical integration method chosen"))
	 , d_oneRotationPerIntegrationPoint(initData(&d_oneRotationPerIntegrationPoint,false,"oneRotationPerIntegrationPoint","if true then computes one rotation per integration point"))
    , lambda(0)
    , mu(0)
    , tetrahedronHandler(NULL)
{
    tetrahedronHandler = new FTCFTetrahedronHandler(this,&tetrahedronInfo);
}

template <class DataTypes> BezierTetrahedralCorotationalFEMForceField<DataTypes>::~BezierTetrahedralCorotationalFEMForceField()
{
    if (tetrahedronHandler) delete tetrahedronHandler;
}

 template <class DataTypes> void BezierTetrahedralCorotationalFEMForceField<DataTypes>::assembleAnisotropicTensors()
 {

 }
template <class DataTypes> void BezierTetrahedralCorotationalFEMForceField<DataTypes>::init()
{
    //	serr << "initializing BezierTetrahedralCorotationalFEMForceField" << sendl;
    this->Inherited::init();

    _topology = this->getContext()->getMeshTopology();
	this->getContext()->get(bezierTetraGeo);

    if ((_topology->getNbTetrahedra()==0) || (!bezierTetraGeo))
    {
        serr << "ERROR(BezierTetrahedralCorotationalFEMForceField): object must have a Tetrahedral Set Topology and a BezierTetrahedronSetGeometryAlgorithms component "<<sendl;
        return;
    }
    updateLameCoefficients();


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
        integrationMethod= AFFINE_ELEMENT_INTEGRATION;
    else if (d_integrationMethod.getValue() == "numerical") 
        integrationMethod= NUMERICAL_INTEGRATION;
    else if (d_integrationMethod.getValue() == "standard") 
        integrationMethod= STANDARD_INTEGRATION;
    else
    {
        serr << "cannot recognize method "<< d_integrationMethod.getValue() << ". Must be either \"analytical\" or \"numerical\"  or \"standard\"" << sendl;
    }
    helper::vector<TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());
    tetrahedronInf.resize(_topology->getNbTetrahedra());

    if (_initialPoints.size() == 0)
    {
        // get restPosition
        const VecCoord& p = this->mstate->read(core::ConstVecCoordId::restPosition())->getValue();
        _initialPoints=p;
    }

    size_t i;


	// precompute the coefficients for handling affine elements
	topology::TetrahedronIndexVector tbi1,tbi2;
	affineStiffnessCoefficientArray.clear();
	
	topology::HighOrderDegreeType degree=bezierTetraGeo->getTopologyContainer()->getDegree();
	size_t nbControlPoints=(degree+1)*(degree+2)*(degree+3)/6;
	sofa::helper::vector<topology::TetrahedronIndexVector> tbiArray;
	tbiArray=bezierTetraGeo->getTopologyContainer()->getTetrahedronIndexArray();

	computeElasticityTensor();
	if (integrationMethod== AFFINE_ELEMENT_INTEGRATION) 
	{
		std::vector<Real> coeffArray(6);
		Mat4x4 coeffMatrix;
		size_t j,k,l,m;
		for (j=0;j<nbControlPoints;j++) {
			tbi1=tbiArray[j];
			for (k=j+1;k<nbControlPoints;k++) {
				tbi2=tbiArray[k];
				coeffMatrix=bezierTetraGeo->getAffineStiffnessCoefficientMatrix(tbi1,tbi2);
				// substract the diagonal terms such that only edge stiffness are used
				for(l=0; l<4; ++l){
					for(m=0; m<4; ++m){
						if (m!=l) {
							coeffMatrix[l][m]-=0.5*(coeffMatrix[l][l]+coeffMatrix[m][m]);
						}
					}
				}
				affineStiffnessCoefficientArray.push_back(coeffMatrix);
			}
		}
	}
	if (integrationMethod== NUMERICAL_INTEGRATION) 
	{
		numericalIntegrationStiffnessDataArray.clear();
		/// get value of integration points0
		topology::NumericalIntegrationDescriptor<Real,4> &nid=bezierTetraGeo->getTetrahedronNumericalIntegrationDescriptor();
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

			std::vector<Vec4> shapeFunctionDerivativeArray;
			for(j=0;j<tbiArray.size();++j) {
				Vec4 deriv=bezierTetraGeo->computeShapeFunctionDerivatives(tbiArray[j],qp.first);
				shapeFunctionDerivativeArray.push_back(deriv);
				Deriv der(deriv[0]-deriv[3],deriv[1]-deriv[3],deriv[2]-deriv[3]);
				nimd.coefficientArray.push_back(der);
			}
			for(j=0;j<tbiArray.size();++j) {
				for(k=j+1;k<tbiArray.size();++k) {
					coeffMatrix=dyad(shapeFunctionDerivativeArray[j],shapeFunctionDerivativeArray[k])*6*weight;
					for(l=0; l<4; ++l){
						for(m=0; m<4; ++m){
							if (l!=m) {
								coeffMatrix[l][m]-=0.5*(coeffMatrix[l][l]+coeffMatrix[m][m]);
							}
						}
					}
					nimd.weightArray.push_back(coeffMatrix);
				}
			}

			numericalIntegrationStiffnessDataArray.push_back(nimd);
		}
	}
	if (integrationMethod== STANDARD_INTEGRATION) 
	{
		numericalIntegrationStiffnessDataArray.clear();
		/// get value of integration points0
		topology::NumericalIntegrationDescriptor<Real,4> &nid=bezierTetraGeo->getTetrahedronNumericalIntegrationDescriptor();
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
				Vec4 deriv=bezierTetraGeo->computeShapeFunctionDerivatives(tbiArray[j],qp.first);
				Deriv der(deriv[0]-deriv[3],deriv[1]-deriv[3],deriv[2]-deriv[3]);
				nimd.coefficientArray.push_back(der);
			}
			numericalIntegrationStiffnessDataArray.push_back(nimd);
		}



	}

	    /// initialize the data structure associated with each tetrahedron
    for (i=0; i<_topology->getNbTetrahedra(); ++i)
    {
        tetrahedronHandler->applyCreateFunction(i,tetrahedronInf[i],_topology->getTetrahedron(i),
                (const helper::vector< unsigned int > )0,
                (const helper::vector< double >)0);
    }
    /// set the call back function upon creation of a tetrahedron
    tetrahedronInfo.createTopologicalEngine(_topology,tetrahedronHandler);
    tetrahedronInfo.registerTopologicalData();
    tetrahedronInfo.endEdit();

	updateTopologyInfo=true;

}


template <class DataTypes>
void BezierTetrahedralCorotationalFEMForceField<DataTypes>::updateTopologyInformation()
{
    int i;
    unsigned int j;

    int nbTetrahedra=_topology->getNbTetrahedra();

    TetrahedronRestInformation *tetinfo;

    helper::vector<typename BezierTetrahedralCorotationalFEMForceField<DataTypes>::TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());
   


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
void BezierTetrahedralCorotationalFEMForceField<DataTypes>::computeElasticityTensor() 											
{
	const helper::vector<Real> & anisotropyParameter=d_anisotropyParameter.getValue();
	const helper::vector<Coord> & anisotropyDirection=d_anisotropyDirection.getValue();
	if (elasticitySymmetry==ISOTROPIC) {
			// elasticity tensor in isotropic case
		Real lambda=getLambda();
		Real mu=getMu();
		elasticityTensor(0,0)=2*mu+lambda;elasticityTensor(1,1)=2*mu+lambda;elasticityTensor(2,2)=2*mu+lambda;
		elasticityTensor(0,1)=lambda;elasticityTensor(0,2)=lambda;elasticityTensor(1,2)=lambda;
		elasticityTensor(1,0)=lambda;elasticityTensor(2,0)=lambda;elasticityTensor(2,1)=lambda;
		elasticityTensor(3,3)=mu;elasticityTensor(4,4)=mu;elasticityTensor(5,5)=mu;
	} else if (elasticitySymmetry==CUBIC) {
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
		Mat3x3 Nd;
		Nd.identity();
		Nd/= sqrt(3.0f);
		Mat3x3 Ne=(2*dyad(n,n)-dyad(v1,v1)-dyad(v2,v2))/sqrt(6.0f);
		Mat3x3 Np=(dyad(v1,v1)-dyad(v2,v2))/sqrt(2.0f);
		Mat3x3 Ns1=(dyad(v2,n)+dyad(n,v2))/sqrt(2.0f);
		Mat3x3 Ns2=(dyad(v1,n)+dyad(n,v1))/sqrt(2.0f);
		Mat3x3 Ns3=(dyad(v1,v2)+dyad(v2,v1))/sqrt(2.0f);
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
		
		Mat3x3 tmp;
		Real val1=(0.5*(1+salpha)+sqrt(2.0)*calpha/4.0f); 
		Real val2=(0.5*(1-salpha)+sqrt(2.0)*calpha/2.0f); 
		Mat3x3 Nh1=val1*(dyad(v1,v1)+dyad(v2,v2))+val2*dyad(n,n);
		Nh1/=sqrt(2*val1*val1+val2*val2);
		val1=(0.5*(1-salpha)-sqrt(2.0)*calpha/4.0f); 
		val2=(0.5*(1+salpha)-sqrt(2.0)*calpha/2.0f); 
		Mat3x3 Nh2=val1*(dyad(v1,v1)+dyad(v2,v2))+val2*dyad(n,n);
		Nh2/=sqrt(2*val1*val1+val2*val2);

		Mat3x3 Np=(dyad(v1,v1)-dyad(v2,v2))/sqrt(2.0f);
		Mat3x3 Ns1=(dyad(v2,n)+dyad(n,v2))/sqrt(2.0f);
		Mat3x3 Ns2=(dyad(v1,n)+dyad(n,v1))/sqrt(2.0f);
		Mat3x3 Ns3=(dyad(v1,v2)+dyad(v2,v1))/sqrt(2.0f);
	
		
		
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
	}
}

template<class DataTypes>
void BezierTetrahedralCorotationalFEMForceField<DataTypes>::computeTetrahedronStiffnessEdgeMatrix(const Coord point[4],
																								  Mat3x3 edgeStiffness[6])
{
	Coord shapeVector[4];
	/// compute 6 times the rest volume
	Real volume=dot(cross(point[1]-point[0],point[2]-point[0]),point[0]-point[3]);
	/// store the rest volume
//	my_tinfo.restVolume=volume/6;

	size_t j,k,l,m,n;
	// store shape vectors at the rest configuration
	for(j=0; j<4; ++j)
	{
		if ((j%2)==0)
			shapeVector[j]=cross(point[(j+2)%4] - point[(j+1)%4],point[(j+3)%4] - point[(j+1)%4])/volume;
		else
			shapeVector[j]= -cross(point[(j+2)%4] - point[(j+1)%4],point[(j+3)%4] - point[(j+1)%4])/volume;
		
	}
	if (elasticitySymmetry==ISOTROPIC) {
		Real mu=getMu()*fabs(volume)/6;
		Real lambda=getLambda()*fabs(volume)/6;
		Real val;

		/// compute the edge stiffness of the linear elastic material
		for(j=0; j<6; ++j)
		{
			k=edgesInTetrahedronArray[j][0];
			l=edgesInTetrahedronArray[j][1];
			// the linear stiffness matrix using shape vectors and Lame coefficients
			val=mu*dot(shapeVector[l],shapeVector[k]);
			for(m=0; m<3; ++m)
			{
				for(n=0; n<3; ++n)
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
		for(j=0; j<6; ++j)
		{
			k=edgesInTetrahedronArray[j][0];
			l=edgesInTetrahedronArray[j][1];
			// the linear stiffness matrix using shape vectors and Lame coefficients
			Mat3x3 tmp=dyad(shapeVector[l],shapeVector[k]);
			for(i=0;i<anisotropyScalarArray.size();++i) {
				edgeStiffness[j]+=anisotropyScalarArray[i]*anisotropyMatrixArray[j]*tmp*anisotropyMatrixArray[j];
			}
		}
	}
}
template<class DataTypes>
void BezierTetrahedralCorotationalFEMForceField<DataTypes>::computeQRRotation( Mat3x3 &r, const Coord *dp)
{
    // first vector on first edge
    // second vector in the plane of the two first edges
    // third vector orthogonal to first and second

    Coord edgex = dp[0];
    edgex.normalize();

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
    r[2][2] = edgez[2];
}

template <class DataTypes>
const helper::vector<typename BezierTetrahedralCorotationalFEMForceField<DataTypes>::Mat3x3> & 
	BezierTetrahedralCorotationalFEMForceField<DataTypes>::getStiffnessArray(
	const size_t i,
	typename BezierTetrahedralCorotationalFEMForceField<DataTypes>::TetrahedronRestInformation *restTetra)
{
	return(restTetra->stiffnessVector);
}
template <class DataTypes>
const  helper::vector<typename BezierTetrahedralCorotationalFEMForceField<DataTypes>::Mat3x3> & 
	BezierTetrahedralCorotationalFEMForceField<DataTypes>::getRotatedStiffnessArray(
	const size_t i,
	const typename BezierTetrahedralCorotationalFEMForceField<DataTypes>::TetrahedronRestInformation *restTetra)
{

	if (decompositionMethod==LINEAR_ELASTIC)
	{
		return(restTetra->stiffnessVector);
	} else
		return(restTetra->rotatedStiffnessVector);
}
template <class DataTypes>
void BezierTetrahedralCorotationalFEMForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* mparams /* PARAMS FIRST */, 
																	 DataVecDeriv &  dataF, const DataVecCoord &  dataX , const DataVecDeriv & /*dataV*/ )
{

    VecDeriv& f        = *(dataF.beginEdit());
    const VecCoord& x  =   dataX.getValue()  ;
    const VecCoord& x0= this->mstate->read(core::ConstVecCoordId::restPosition())->getValue();
	sofa::helper::vector<Coord> dp,force;
    size_t i,j,k,l,v0,v1,rank,p;
    size_t nbTetrahedra=_topology->getNbTetrahedra();
    
	HighOrderDegreeType degree=bezierTetraGeo->getTopologyContainer()->getDegree();
	size_t nbControlPoints=(degree+1)*(degree+2)*(degree+3)/6;
	HighOrderTetrahedronSetTopologyContainer::VecPointID indexArray;


    if (updateTopologyInfo)
    {
        updateTopologyInformation();
    }
    helper::vector<TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());
    TetrahedronRestInformation *tetinfo;
	
	dp.resize(nbControlPoints);
	force.resize(nbControlPoints);
	
	Coord dpos,sv;
		
	for(i=0; i<nbTetrahedra; i++ )
	{
		tetinfo=&tetrahedronInf[i];
		const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=
			bezierTetraGeo->getTopologyContainer()->getGlobalIndexArrayOfControlPoints(i);
		
		nbControlPoints=indexArray.size();
	
		if (d_oneRotationPerIntegrationPoint.getValue()) {
			
			Mat3x3 S,R;
			size_t j,k,l,m,n;
//			helper::vector<Mat3x3> stiffnessArray(nbControlPoints*(nbControlPoints-1)/2);
			helper::vector<Mat3x3> &stiffnessArray=tetinfo->rotatedStiffnessVector;
			std::fill(stiffnessArray.begin(),stiffnessArray.end(),Mat3x3());
			assert(stiffnessArray.size()==nbControlPoints*(nbControlPoints-1)/2);
			// loop through the integration points
			for (l=0;l<numericalIntegrationStiffnessDataArray.size();++l) {
				Coord dpp[6],point[4];
				// the barycentric coordinate
				bezierTetraGeo->computeNodalValueDerivatives(i,numericalIntegrationStiffnessDataArray[l].integrationPoint, x,point);

				for(j=0; j<6; ++j){
					m=edgesInTetrahedronArray[j][0];
					n=edgesInTetrahedronArray[j][1];

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
//				Real w=numericalIntegrationStiffnessDataArray[l].integrationWeight;
				for (rank=0,j=0; j<nbControlPoints; ++j) {
					v0 = indexArray[j];
					for ( k=j+1; k<nbControlPoints; ++k,++rank) {

						// compute K[j][k] from the edge stiffness array stored in integrationPointsStiffnessVector
						Mat3x3 stiffness;


						const Mat4x4  & coeffMatrix=numericalIntegrationStiffnessDataArray[l].weightArray[rank];
						/// add edge stiffness
						for(p=0;p<6;++p) {
							m=edgesInTetrahedronArray[p][0];
							n=edgesInTetrahedronArray[p][1];
							const Mat3x3  &edgeStiff =tetinfo->integrationPointsStiffnessVector[l][p];

							stiffness+= edgeStiff*coeffMatrix[m][n]+
								edgeStiff.transposed()*coeffMatrix[n][m];
						}



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
							Mat3x3 mat=R*stiffness*R.transposed();
							stiffnessArray[rank]+=mat;
						}
//				if (rank==0) {
//					std::cerr<< stiffness<< "R="<<R<<" stiffArray="<<stiffnessArray[0]<<std::endl;
//				}
					}

				}
			}
//	tetinfo->rotatedStiffnessVector.clear();
//	tetinfo->rotatedStiffnessVector=stiffnessArray;
			for (j=0; j<nbControlPoints; ++j) {
				f[indexArray[j]]+=R*force[j];
			}


		} else {
			const  helper::vector<Mat3x3> &stiffnessArray=getStiffnessArray(i,tetinfo);
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
				Mat3x3 deformationGradient,S,R;
				Coord dpp[6];
				for (j=0; j<6; ++j)
				{
					dpp[j]=x[tetinfo->v[edgesInTetrahedronArray[j][1]]]-x[tetinfo->v[edgesInTetrahedronArray[j][0]]];
				}
				if (decompositionMethod==POLAR_DECOMPOSITION)
				{
					// compute the deformation gradient
					// deformation gradient = sum of tensor product between vertex position and shape vector
					// optimize by using displacement with first vertex
					sv=tetinfo->shapeVector[1];


					for (k=0; k<3; ++k)
					{
						for (l=0; l<3; ++l)
						{
							deformationGradient[k][l]=dpp[0][k]*sv[l];
						}
					}
					for (j=1; j<3; ++j)
					{
						sv=tetinfo->shapeVector[j+1];
						for (k=0; k<3; ++k)
						{
							for (l=0; l<3; ++l)
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
							Mat3x3 mat=R*stiffnessArray[rank]*tetinfo->rotation;
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
    tetrahedronInfo.endEdit();

    dataF.endEdit();

}


template <class DataTypes>
void BezierTetrahedralCorotationalFEMForceField<DataTypes>::addDForce(const sofa::core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv&   datadF , const DataVecDeriv&   datadX )
{
    VecDeriv& df       = *(datadF.beginEdit());
    const VecCoord& dx =   datadX.getValue()  ;
    Real kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());
    const helper::vector<TetrahedronRestInformation>& tetrahedronInf = tetrahedronInfo.getValue();
	Coord dpos;
    size_t i,j,k,v0,v1,rank;
    size_t nbTetrahedra=_topology->getNbTetrahedra();
    
	HighOrderDegreeType degree=bezierTetraGeo->getTopologyContainer()->getDegree();
	size_t nbControlPoints=(degree+1)*(degree+2)*(degree+3)/6;
	
    const TetrahedronRestInformation *tetinfo;
	

	for(i=0; i<nbTetrahedra; i++ )
	{
		tetinfo=&tetrahedronInf[i];
		const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=
			bezierTetraGeo->getTopologyContainer()->getGlobalIndexArrayOfControlPoints(i);
		const  helper::vector<Mat3x3> &stiffnessArray=getRotatedStiffnessArray(i,tetinfo);
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
		/*
		indexArray.clear();
		/// get the global index of each control point in the tetrahedron
		bezierTetraGeo->getTopologyContainer()->getGlobalIndexArrayOfBezierPointsInTetrahedron(i,indexArray) ;
		if (decompositionMethod==LINEAR_ELASTIC)
		{
			// loop over each entry in the stiffness vector of size nbControlPoints*(nbControlPoints+1)/2
			for (rank=0,j=0; j<nbControlPoints; ++j) {
				v0 = indexArray[j];
				for ( k=j+1; k<nbControlPoints; ++k,++rank) {
					v1 = indexArray[k];
					dpos=dx[v0]-dx[v1];
				
					df[v1]-=tetinfo->stiffnessVector[rank]*dpos*kFactor;
					df[v0]+=tetinfo->stiffnessVector[rank].multTranspose(dpos*kFactor);
				}
			}
		} else {
			Mat3x3 rot=tetinfo->rotation;
			// loop over each entry in the stiffness vector of size nbControlPoints*(nbControlPoints+1)/2
			for (rank=0,j=0; j<nbControlPoints; ++j) {
				v0 = indexArray[j];
				for ( k=j+1; k<nbControlPoints; ++k,++rank) {
					v1 = indexArray[k];
					dpos=dx[v0]-dx[v1];
					df[v1]-=rot.multTranspose(tetinfo->stiffnessVector[rank]*rot*dpos*kFactor);
					df[v0]+=rot.multTranspose(tetinfo->stiffnessVector[rank].multTranspose(rot*dpos*kFactor));
				}
			}

		}
		*/
	}

    datadF.endEdit();
}

template<class DataTypes>
SReal BezierTetrahedralCorotationalFEMForceField<DataTypes>::getPotentialEnergy(const core::MechanicalParams*, const DataVecCoord&) const
{
    serr << "ERROR("<<this->getClassName()<<"): getPotentialEnergy( const MechanicalParams*, const DataVecCoord& ) not implemented." << sendl;
    return 0.0;
}


template<class DataTypes>
void BezierTetrahedralCorotationalFEMForceField<DataTypes>::updateLameCoefficients()
{
    lambda= d_youngModulus.getValue()*d_poissonRatio.getValue()/((1-2*d_poissonRatio.getValue())*(1+d_poissonRatio.getValue()));
    mu = d_youngModulus.getValue()/(2*(1+d_poissonRatio.getValue()));

}


template<class DataTypes>
void BezierTetrahedralCorotationalFEMForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
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

#endif // SOFA_COMPONENT_FORCEFIELD_BezierTetrahedralCorotationalFEMForceField_INL
