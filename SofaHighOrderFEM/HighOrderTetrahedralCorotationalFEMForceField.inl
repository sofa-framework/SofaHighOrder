
#ifndef SOFA_COMPONENT_FORCEFIELD_HIGHORDERTETRAHEDRALCOROTATIONALFEMFORCEFIELD_INL
#define SOFA_COMPONENT_FORCEFIELD_HIGHORDERTETRAHEDRALCOROTATIONALFEMFORCEFIELD_INL

#include "HighOrderTetrahedralCorotationalFEMForceField.h"
#include <sofa/core/visual/VisualParams.h>
#include <fstream> // for reading the file
#include <iostream> //for debugging
#include <sofa/helper/gl/template.h>
#include <SofaBaseTopology/TopologyData.inl>
#include <HighOrderTetrahedronSetGeometryAlgorithms.h>
#include <BezierTetrahedronSetGeometryAlgorithms.h>
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

typedef BaseMeshTopology::Tetra				Tetra;
typedef BaseMeshTopology::EdgesInTetrahedron		EdgesInTetrahedron;

typedef Tetra			        Tetrahedron;
typedef EdgesInTetrahedron		EdgesInTetrahedron;

const unsigned int edgesInTetrahedronArray[6][2] = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}};




template< class DataTypes>
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::weightArrayPointer::allocate(size_t degree) {
	if (degree==2) {
		weightArrayQuadratic[0]=boost::make_shared<Mat45x6>();
		weightArrayQuadratic[1]=boost::make_shared<Mat45x6>();

	} else 	if (degree==3) {
		weightArrayCubic[0]=boost::make_shared<Mat190x6>();
		weightArrayCubic[1]=boost::make_shared<Mat190x6>();
	} else 	if (degree==4) {
		weightArrayQuartic[0]=boost::make_shared<Mat595x6>();
		weightArrayQuartic[1]=boost::make_shared<Mat595x6>();
	} else 	if (degree==5) {
		weightArrayQuintic[0]=boost::make_shared<Mat1540x6>();
		weightArrayQuintic[1]=boost::make_shared<Mat1540x6>();
	}
}


template <int L, class real = float>
Mat<L, L, real> symmetrizeMatrix(Mat<L, L, real> m)  {
    size_t i,j;
    Mat<L, L, real> res=m;
    for (i = 0; i < L; ++i) {
        res[i][i] = 2 * m[i][i];
        for (j = i+1; j < L; ++j) {
            res[i][j] = m[i][j] + m[j][i];
            res[j][i] = res[i][j];
        }
    }
    return res;
}


template< class DataTypes>
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::FTCFTetrahedronHandler::applyCreateFunction(unsigned int tetrahedronIndex,
    TetrahedronRestInformation &my_tinfo,
    const Tetrahedron &tt,
    const sofa::helper::vector<unsigned int> &,
    const sofa::helper::vector<double> &)
{
    if (ff)
    {


        const std::vector< Tetrahedron > &tetrahedronArray = ff->_topology->getTetrahedra();
        HighOrderTetrahedronSetTopologyContainer *container = ff->highOrderTetraGeo->getTopologyContainer();
        HighOrderDegreeType degree = container->getDegree();
        size_t nbControlPoints = (degree + 1)*(degree + 2)*(degree + 3) / 6;
        size_t nbStiffnessEntries = nbControlPoints*(nbControlPoints - 1) / 2;
        if (my_tinfo.stiffnessVector.size() != nbStiffnessEntries) {
            my_tinfo.stiffnessVector.resize(nbStiffnessEntries);
        }
        // set array to zero
        std::fill(my_tinfo.stiffnessVector.begin(), my_tinfo.stiffnessVector.end(), Mat3x3());

        //		const std::vector< Edge> &edgeArray=ff->_topology->getEdges() ;
        size_t i, j, k, l, m, n;

        typename DataTypes::Coord point[4];


        const typename DataTypes::VecCoord &restPosition = ff->mstate->read(core::ConstVecCoordId::restPosition())->getValue();
        // now computed the stiffness for the HighOrder Tetrahedron
        sofa::helper::vector<TetrahedronIndexVector> tbiArray;

        tbiArray = ff->highOrderTetraGeo->getTopologyContainer()->getTetrahedronIndexArray();

        size_t rank = 0;

        if (ff->d_oneRotationPerIntegrationPoint.getValue()) {
            // one rotation per integration point
            Mat6x9 edgeStiffness[2];  // the off-diagonal 3x3 block matrices that makes the 12x12 linear elastic matrix
            sofa::defaulttype::Vec<4, Real> bc;
            //			my_tinfo.integrationPointsStiffnessVector.resize(ff->numericalIntegrationStiffnessDataArray.size());
            my_tinfo.integrationPointsRestEdgeVector.resize(ff->numericalIntegrationStiffnessDataArray.size());
            my_tinfo.integrationPointsRestRotationArray.resize(ff->numericalIntegrationStiffnessDataArray.size());

            // loop through the integration points
            for (i = 0; i < ff->numericalIntegrationStiffnessDataArray.size(); ++i) {

                // the barycentric coordinate
                bc = ff->numericalIntegrationStiffnessDataArray[i].integrationPoint;
                // Compute the local tetrahedron by storing in point the 4 derivatives of the shape functions  
                ff->highOrderTetraGeo->computeNodalValueDerivatives(tetrahedronIndex, bc, restPosition, point);

                // initialize rotation of element
                if (ff->decompositionMethod == QR_DECOMPOSITION) {
                    helper::vector<Coord> restEdgeVector(6);
                    Coord restEdgeVectorTmp[6];

                    for (j = 0; j < 6; ++j) {
                        k = edgesInTetrahedronArray[j][0];
                        l = edgesInTetrahedronArray[j][1];

                        // store the rest edge vector
                        restEdgeVector[j] = point[l] - point[k];
                        restEdgeVectorTmp[j] = restEdgeVector[j];
                    }
                    my_tinfo.integrationPointsRestEdgeVector[i] = restEdgeVector;

                    // compute the rotation matrix of the initial tetrahedron for the QR decomposition
                    computeQRRotation(my_tinfo.integrationPointsRestRotationArray[i],
                        restEdgeVectorTmp);
                }
                else 	if (ff->decompositionMethod == POLAR_DECOMPOSITION_MODIFIED) {
                    Mat3x3 Transformation;
                    Transformation[0] = point[1] - point[0];
                    Transformation[1] = point[2] - point[0];
                    Transformation[2] = point[3] - point[0];
                    helper::Decompose<Real>::polarDecomposition(Transformation, my_tinfo.integrationPointsRestRotationArray[i]);
                }


                // compute the edge stiffness associated with that local tetrahedron
                ff->computeTetrahedronStiffnessEdgeMatrix(point, edgeStiffness);


                my_tinfo.integrationPointsStiffnessVector.push_back(edgeStiffness[0]);
                my_tinfo.integrationPointsStiffnessVector.push_back(edgeStiffness[1]);

                my_tinfo.rotatedStiffnessVector.resize(nbControlPoints*(nbControlPoints - 1) / 2);
            }

        }
        else {

            ///describe the indices of the 4 tetrahedron vertices
            const Tetrahedron &t = tetrahedronArray[tetrahedronIndex];
            //    BaseMeshTopology::EdgesInTetrahedron te=ff->_topology->getEdgesInTetrahedron(tetrahedronIndex);


            // store the point position
            for (j = 0; j < 4; ++j)
                point[j] = (restPosition)[t[j]];

            if (ff->decompositionMethod == QR_DECOMPOSITION) {
                for (j = 0; j < 6; ++j) {
                    k = edgesInTetrahedronArray[j][0];
                    l = edgesInTetrahedronArray[j][1];

                    // store the rest edge vector
                    my_tinfo.restEdgeVector[j] = point[l] - point[k];
                }
                // compute the rotation matrix of the initial tetrahedron for the QR decomposition
                computeQRRotation(my_tinfo.restRotation, my_tinfo.restEdgeVector);
            }
            else 	if (ff->decompositionMethod == POLAR_DECOMPOSITION_MODIFIED) {
                Mat3x3 Transformation;
                Transformation[0] = point[1] - point[0];
                Transformation[1] = point[2] - point[0];
                Transformation[2] = point[3] - point[0];
                helper::Decompose<Real>::polarDecomposition(Transformation, my_tinfo.restRotation);
            }

            if ((ff->integrationMethod == HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::AFFINE_ELEMENT_INTEGRATION) ||
                ((ff->d_forceAffineAssemblyForAffineElements.getValue()) && ((ff->highOrderTetraGeo->isBezierTetrahedronAffine(tetrahedronIndex, restPosition))))) {
                Mat6x9 edgeStiffnessVectorized[2];
                ff->computeTetrahedronStiffnessEdgeMatrix(point, edgeStiffnessVectorized);
                helper::system::thread::ctime_t startUpdateMat = helper::system::thread::CTime::getTime();
                if (degree == 1) {
                    for (rank = 0; rank < nbStiffnessEntries; rank++) {
                        my_tinfo.stiffnessVector[rank] += Mat3x3((const Real *)&edgeStiffnessVectorized[0][rank][0]);
                    }
                }
                else if (degree == 2) {
                    Mat45x9 res = (*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayQuadratic[0]))*edgeStiffnessVectorized[0] +
                        (*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayQuadratic[1]))*edgeStiffnessVectorized[1];
                    for (rank = 0; rank < nbStiffnessEntries; rank++) {
                        my_tinfo.stiffnessVector[rank] += Mat3x3((const Real *)&res[rank][0]);
                    }
                }
                else if (degree == 3) {
                    Mat190x9 res = (*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayCubic[0]))*edgeStiffnessVectorized[0] +
                        (*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayCubic[1]))*edgeStiffnessVectorized[1];
                    for (rank = 0; rank < nbStiffnessEntries; rank++) {
                        my_tinfo.stiffnessVector[rank] += Mat3x3((const Real *)&res[rank][0]);
                    }

                }
                else if (degree == 4) {
                    Mat595x9 res = (*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayQuartic[0]))*edgeStiffnessVectorized[0] +
                        (*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayQuartic[1]))*edgeStiffnessVectorized[1];
                    for (rank = 0; rank < nbStiffnessEntries; rank++) {
                        my_tinfo.stiffnessVector[rank] += Mat3x3((const Real *)&res[rank][0]);
                    }
                }
                else if (degree == 5) {
                    Mat1540x9 res = (*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayQuintic[0]))*edgeStiffnessVectorized[0] +
                        (*(ff->affineStiffnessCoefficientPreStoredArray.weightArrayQuintic[1]))*edgeStiffnessVectorized[1];
                    for (rank = 0; rank < nbStiffnessEntries; rank++) {
                        my_tinfo.stiffnessVector[rank] += Mat3x3((const Real *)&res[rank][0]);
                    }
                }
                else {
                    for (rank = 0; rank < nbStiffnessEntries; rank++) {
                        const Vec6  & coeffVec1 = ff->affineStiffnessCoefficientArray[2 * rank];
                        const Vec6  & coeffVec2 = ff->affineStiffnessCoefficientArray[2 * rank + 1];

                        //			Vec9 res=edgeStiffness[0]*coeffVec1+edgeStiffness[1]*coeffVec2;
                        Vec9 res = edgeStiffnessVectorized[0].multTranspose(coeffVec1) + edgeStiffnessVectorized[1].multTranspose(coeffVec2);
                        my_tinfo.stiffnessVector[rank] += Mat3x3((const Real *)&res[0]);
                    }
                }

                if (ff->f_printLog.getValue()) {
                    helper::system::thread::ctime_t endUpdateMat = helper::system::thread::CTime::getTime();
                    ff->totalUpdateMat += endUpdateMat - startUpdateMat;
                }

            }
            else if (ff->integrationMethod == HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::STANDARD_INTEGRATION) {
                sofa::defaulttype::Vec<4, Real> bc;

                // loop through the integration points
                for (i = 0; i < ff->numericalIntegrationStiffnessDataArray.size(); ++i) {

                    // the barycentric coordinate
                    bc = ff->numericalIntegrationStiffnessDataArray[i].integrationPoint;
                    // Compute the local tetrahedron by storing in point the 4 derivatives of the shape functions  
                    ff->highOrderTetraGeo->computeNodalValueDerivatives(tetrahedronIndex, bc, restPosition, point);

                    Mat3x3 Jacobian, inverseJacobian;
                    for (j = 0; j < 3; ++j) {
                        for (k = 0; k < 3; ++k) {
                            Jacobian[j][k] = point[j][k] - point[3][k];
                        }
                    }
                    invertMatrix(inverseJacobian, Jacobian);
                    Real jac = fabs(determinant(Jacobian))*ff->numericalIntegrationStiffnessDataArray[i].integrationWeight;
                    helper::system::thread::ctime_t startUpdateMat = helper::system::thread::CTime::getTime();
                    helper::vector<Mat6x3> SDArray;
                    for (j = 0; j < nbControlPoints; j++) {
                        Coord sv = inverseJacobian*ff->numericalIntegrationStiffnessDataArray[i].coefficientArray[j];
                        Mat6x3 strainDisplacement;
                        strainDisplacement[0][0] = sv[0]; strainDisplacement[1][1] = sv[1]; strainDisplacement[2][2] = sv[2];
                        strainDisplacement[3][1] = sv[0]; strainDisplacement[3][0] = sv[1]; strainDisplacement[4][1] = sv[2];
                        strainDisplacement[5][2] = sv[0]; strainDisplacement[4][2] = sv[1]; strainDisplacement[5][0] = sv[2];
                        SDArray.push_back(strainDisplacement);
                    }

                    for (rank = 0, j = 0; j < nbControlPoints; j++) {

                        for (k = j + 1; k < nbControlPoints; k++, rank++) {
                            Mat3x3 edgeStiffness = (SDArray[j].transposed()*((ff->elasticityTensor)*SDArray[k]));
                            edgeStiffness *= jac;
                            my_tinfo.stiffnessVector[rank] += edgeStiffness.transposed();
                        }
                    }
                    if (ff->f_printLog.getValue()) {
                        helper::system::thread::ctime_t endUpdateMat = helper::system::thread::CTime::getTime();
                        ff->totalUpdateMat += endUpdateMat - startUpdateMat;
                    }
                }

            }
            else if (ff->integrationMethod == HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::NUMERICAL_INTEGRATION) {
                sofa::defaulttype::Vec<4, Real> bc;

                size_t p, q;
                Mat6x9 edgeStiffnessVectorized[2];
                // loop through the integration points
                for (i = 0; i < ff->numericalIntegrationStiffnessDataArray.size(); ++i) {
                    // copy weight array locally to speed-up computation
                    //		std::vector<Mat4x4>  weightArray=ff->numericalIntegrationStiffnessDataArray[i].weightArray;
                    // the barycentric coordinate
                    bc = ff->numericalIntegrationStiffnessDataArray[i].integrationPoint;
                    // Compute the local tetrahedron by storing in point the 4 derivatives of the shape functions  
                    ff->highOrderTetraGeo->computeNodalValueDerivatives(tetrahedronIndex, bc, restPosition, point);
                    // compute the edge stiffness associated with that local tetrahedron
                    ff->computeTetrahedronStiffnessEdgeMatrix(point, edgeStiffnessVectorized);
                    helper::system::thread::ctime_t startUpdateMat = helper::system::thread::CTime::getTime();
                    // compute the stiffness matrix for all pairs of control points 

                    if (degree == 2) {
                        Mat45x9 res = (*(ff->numericalIntegrationStiffnessDataArray[i].arrayPointer.weightArrayQuadratic[0]))*edgeStiffnessVectorized[0] +
                            (*(ff->numericalIntegrationStiffnessDataArray[i].arrayPointer.weightArrayQuadratic[1]))*edgeStiffnessVectorized[1];
                        for (rank = 0; rank < nbStiffnessEntries; rank++) {
                            my_tinfo.stiffnessVector[rank] += Mat3x3((const Real *)&res[rank][0]);
                        }
                    }
                    else if (degree == 3) {
                        Mat190x9 res = (*(ff->numericalIntegrationStiffnessDataArray[i].arrayPointer.weightArrayCubic[0]))*edgeStiffnessVectorized[0] +
                            (*(ff->numericalIntegrationStiffnessDataArray[i].arrayPointer.weightArrayCubic[1]))*edgeStiffnessVectorized[1];
                        for (rank = 0; rank < nbStiffnessEntries; rank++) {
                            my_tinfo.stiffnessVector[rank] += Mat3x3((const Real *)&res[rank][0]);
                        }
                    }
                    else if (degree == 4) {
                        Mat595x9 res = (*(ff->numericalIntegrationStiffnessDataArray[i].arrayPointer.weightArrayQuartic[0]))*edgeStiffnessVectorized[0] +
                            (*(ff->numericalIntegrationStiffnessDataArray[i].arrayPointer.weightArrayQuartic[1]))*edgeStiffnessVectorized[1];
                        for (rank = 0; rank < nbStiffnessEntries; rank++) {
                            my_tinfo.stiffnessVector[rank] += Mat3x3((const Real *)&res[rank][0]);
                        }
                    }
                    else if (degree == 5) {
                        Mat1540x9 res = (*(ff->numericalIntegrationStiffnessDataArray[i].arrayPointer.weightArrayQuintic[0]))*edgeStiffnessVectorized[0] +
                            (*(ff->numericalIntegrationStiffnessDataArray[i].arrayPointer.weightArrayQuintic[1]))*edgeStiffnessVectorized[1];
                        for (rank = 0; rank < nbStiffnessEntries; rank++) {
                            my_tinfo.stiffnessVector[rank] += Mat3x3((const Real *)&res[rank][0]);
                        }
                    }
                    else {
                        const std::vector<Vec6> & weightArray = ff->numericalIntegrationStiffnessDataArray[i].weightArray;
                        for (rank = 0; rank < nbStiffnessEntries; rank++) {
                            const Vec6  & coeffVec1 = weightArray[2 * rank];
                            const Vec6  & coeffVec2 = weightArray[2 * rank + 1];

                            //			Vec9 res=edgeStiffness[0]*coeffVec1+edgeStiffness[1]*coeffVec2;
                            Vec9 res = edgeStiffnessVectorized[0].multTranspose(coeffVec1) + edgeStiffnessVectorized[1].multTranspose(coeffVec2);
                            my_tinfo.stiffnessVector[rank] += Mat3x3((const Real *)&res[0]);
                        }
                    }

                    if (ff->f_printLog.getValue()) {
                        helper::system::thread::ctime_t endUpdateMat = helper::system::thread::CTime::getTime();
                        ff->totalUpdateMat += endUpdateMat - startUpdateMat;
                    }

                }
            }
            else if (ff->integrationMethod == HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::BEZIER_NUMERICAL_INTEGRATION) {
                helper::system::thread::ctime_t startUpdateMat = helper::system::thread::CTime::getTime();
                sofa::defaulttype::Vec<4, Real> bc;
                Mat6x9 edgeStiffnessVectorized[2];
                std::vector<  sofa::defaulttype::Mat<6, 9, Real> > reducedStiffness;
                
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
                    ff->computeTetrahedronStiffnessEdgeMatrix(point, edgeStiffnessVectorized);
                    for (j = 0; j < numberReducedEntries; ++j) {
                        reducedStiffness[j] += ff->numericalIntegrationStiffnessDataArray[i].weightBezierArray[j] * edgeStiffnessVectorized[0];
                      }
                }
                size_t r;
                for (i = 0; i < nbStiffnessEntries; ++i) {
                    sofa::defaulttype::Vec<16, int>  &mapping = ff->bezierMappingArray[i];
                    for (rank = 0, j = 0; j < 4; ++j) {
                        for (k = j + 1; k < 4; ++k, ++rank) {
                            r = j * 4 + k;
                            if (mapping[r] >= 0)
                                my_tinfo.stiffnessVector[i] += ff->bezierCoefficientArray[i][r] * Mat3x3(&reducedStiffness[(size_t)mapping[r]][rank][0]);
                            r = k * 4 + j;
                            if (mapping[r] >= 0)
                                my_tinfo.stiffnessVector[i] += ff->bezierCoefficientArray[i][r] * Mat3x3(&reducedStiffness[(size_t)mapping[r]][rank][0]).transposed();
                            r = j * 4 + j;
                            if (mapping[r] >= 0)
                                my_tinfo.stiffnessVector[i] -= 0.5*ff->bezierCoefficientArray[i][r] * symmetrizeMatrix<3, Real>(Mat3x3(&reducedStiffness[(size_t)mapping[r]][rank][0]));
                            r = k * 4 + k;
                            if (mapping[r] >= 0)
                                my_tinfo.stiffnessVector[i] -= 0.5*ff->bezierCoefficientArray[i][r] * symmetrizeMatrix<3, Real>(Mat3x3(&reducedStiffness[(size_t)mapping[r]][rank][0]));
                        }
                    }
                }
                if (ff->f_printLog.getValue()) {
                    helper::system::thread::ctime_t endUpdateMat = helper::system::thread::CTime::getTime();
                    ff->totalUpdateMat += endUpdateMat - startUpdateMat;
                }
            }
            else if (ff->integrationMethod == HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::NUMERICAL_INTEGRATION_2) {
                sofa::defaulttype::Vec<4, Real> bc;
                Mat3x3 edgeStiffness3x3[6];
                size_t p, q;

                // loop through the integration points
                for (i = 0; i < ff->numericalIntegrationStiffnessDataArray.size(); ++i) {
                    // copy weight array locally to speed-up computation
                    //		std::vector<Mat4x4>  weightArray=ff->numericalIntegrationStiffnessDataArray[i].weightArray;
                    // the barycentric coordinate
                    bc = ff->numericalIntegrationStiffnessDataArray[i].integrationPoint;
                    // Compute the local tetrahedron by storing in point the 4 derivatives of the shape functions  
                    ff->highOrderTetraGeo->computeNodalValueDerivatives(tetrahedronIndex, bc, restPosition, point);
                    // compute the edge stiffness associated with that local tetrahedron
                    ff->computeTetrahedronStiffnessEdgeMatrix(point, edgeStiffness3x3);
                    helper::system::thread::ctime_t startUpdateMat = helper::system::thread::CTime::getTime();
                    // compute the stiffness matrix for all pairs of control points 

                    const std::vector<Vec6> & weightArray = ff->numericalIntegrationStiffnessDataArray[i].weightArray;



                    for (rank = 0, j = 0; j < nbControlPoints; ++j) {
                        for (k = j + 1; k < nbControlPoints; k++, rank++) {
                            //	Mat3x3 stiff;
                            //	const Vec6  & coeffVec1=ff->numericalIntegrationStiffnessDataArray[i].weightArray[2*rank];
                            //	const Vec6  & coeffVec2=ff->numericalIntegrationStiffnessDataArray[i].weightArray[2*rank+1];
                            Mat4x4  coeffMatrix = ff->numericalIntegrationStiffnessDataArray[i].weightArray4x4[rank];
                            /// add edge stiffness
                            //	Vec9 res=edgeStiffness[0]*coeffVec1+edgeStiffness[1]*coeffVec2;
                            //		my_tinfo.stiffnessVector[rank]+= Mat3x3((const Real *) &res[0]);

                            for (l = 0; l < 6; ++l) {
                                m = edgesInTetrahedronArray[l][0];
                                n = edgesInTetrahedronArray[l][1];
                                /*
                                Real a=coeffMatrix[m][n];
                                Real b=coeffMatrix[n][m];
                                for (p=0;p<3;++p) {
                                my_tinfo.stiffnessVector[rank][p][p]+=edgeStiffness[l][p][p]*(a+b);
                                for (q=p+1;q<3;++q) {
                                my_tinfo.stiffnessVector[rank][p][q]+=edgeStiffness[l][p][q]*a+edgeStiffness[l][q][p]*b;
                                my_tinfo.stiffnessVector[rank][q][p]+=edgeStiffness[l][q][p]*a+edgeStiffness[l][p][q]*b;
                                }
                                }  */
                                //	if (coeffMatrix[0][l]!=0)
                                //		my_tinfo.stiffnessVector[rank]+=edgeStiffness[l]*coeffMatrix[0][l];
                                //	if (coeffMatrix[1][l]!=0)
                                //		my_tinfo.stiffnessVector[rank]+=edgeStiffness[l].transposed()*coeffMatrix[1][l];

                                if (coeffMatrix[m][n] != 0) {
                                    my_tinfo.stiffnessVector[rank] += edgeStiffness3x3[l] * coeffMatrix[m][n];
                                }
                                if (coeffMatrix[n][m] != 0) {
                                    my_tinfo.stiffnessVector[rank] += edgeStiffness3x3[l].transposed()*coeffMatrix[n][m];
                                }

                                //		my_tinfo.stiffnessVector[rank]+= edgeStiffness[l].transposed()*coeffMatrix[n][m] +edgeStiffness[l].transposed()*coeffMatrix[n][m];
                                //	}
                                //	my_tinfo.stiffnessVector[rank]+=stiff;
                            }
                        }
                    }
                    if (ff->f_printLog.getValue()) {
                        helper::system::thread::ctime_t endUpdateMat = helper::system::thread::CTime::getTime();
                        ff->totalUpdateMat += endUpdateMat - startUpdateMat;
                    }
                }


                //std::cerr<<" total update mat="<<((totalAssembly)/(Real)helper::system::thread::CTime::getRefTicksPerSec()) << std::endl;

            }
#ifdef _DEBUG
            if (ff->f_printLog.getValue()) {
                std::cerr << " lambda=" << ff->getLambda() << " mu=" << ff->getMu() << std::endl;


                for (rank = 0, l = 0; l < tbiArray.size(); ++l)
                {
                    for (m = l + 1; m < tbiArray.size(); ++m, ++rank)
                    {
                        std::cerr << "Stiffness entry [" << (unsigned int)tbiArray[l][0] << " " << (unsigned int)tbiArray[l][1] << " " << (unsigned int)tbiArray[l][2] << " " << (unsigned int)tbiArray[l][3] << "][" <<
                            (unsigned int)tbiArray[m][0] << " " << (unsigned int)tbiArray[m][1] << " " << (unsigned int)tbiArray[m][2] << " " << (unsigned int)tbiArray[m][3] << "]=" << my_tinfo.stiffnessVector[rank] << std::endl;

                    }
                }
            }
#endif
            }
        }

    }


template <class DataTypes>
HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::HighOrderTetrahedralCorotationalFEMForceField()
        : tetrahedronInfo(initData(&tetrahedronInfo, "tetrahedronInfo", "Internal tetrahedron data"))
        , _initialPoints(0)
        , updateMatrix(true)
        , d_method(initData(&d_method, std::string("linear"), "method", "method for rotation computation :\"qr\" (by QR) or \"polar\" or \"polar2\" or \"none\" (Linear elastic)"))
        , d_poissonRatio(initData(&d_poissonRatio, (Real)0.3, "poissonRatio", "Poisson ratio in Hooke's law"))
        , d_youngModulus(initData(&d_youngModulus, (Real)1000., "youngModulus", "Young modulus in Hooke's law"))
        , d_anisotropy(initData(&d_anisotropy, std::string("isotropic"), "elasticitySymmetry", "the type of anisotropy for the elasticity tensor :\"isotropic\"  or \"transverseIsotropic\" or \"orthotropic\" or \"cubic\" "))
        , d_anisotropyParameter(initData(&d_anisotropyParameter, "anisotropyParameters", "the elastic parameters for anisotropic materials "))
        , d_anisotropyDirection(initData(&d_anisotropyDirection, "anisotropyDirections", "the directions of anisotropy"))
        , numericalIntegrationOrder(initData(&numericalIntegrationOrder, (size_t)2, "integrationOrder", "The order of integration for numerical integration"))
        , d_integrationMethod(initData(&d_integrationMethod, std::string("analytical"), "integrationMethod", "\"analytical\" if closed form expression for affine element, \"numerical\" if numerical integration is chosen,  \"standard\" if standard integration is chosen"))
        , numericalIntegrationMethod(initData(&numericalIntegrationMethod, std::string("Tetrahedron Gauss"), "numericalIntegrationMethod", "The type of numerical integration method chosen"))
        , d_oneRotationPerIntegrationPoint(initData(&d_oneRotationPerIntegrationPoint, false, "oneRotationPerIntegrationPoint", "if true then computes one rotation per integration point"))
        , d_assemblyTime(initData(&d_assemblyTime, (Real)0, "assemblyTime", "the time spent in assembling the stiffness matrix. Only updated if printLog is set to true"))
        , d_forceAffineAssemblyForAffineElements(initData(&d_forceAffineAssemblyForAffineElements, true, "forceAffineAssemblyForAffineElements", "if true affine tetrahedra are always assembled with the closed form formula, Otherwise use the method defined in integrationMethod"))
        , lambda(0)
        , mu(0)
        , tetrahedronHandler(NULL)
    {
        tetrahedronHandler = new FTCFTetrahedronHandler(this, &tetrahedronInfo);
    }


template <class DataTypes>
HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::~HighOrderTetrahedralCorotationalFEMForceField()
    {
        if (tetrahedronHandler) delete tetrahedronHandler;
    }


template <class DataTypes>
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::assembleAnisotropicTensors()
 {

 }


template <class DataTypes>
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::init()
{
    //	serr << "initializing HighOrderTetrahedralCorotationalFEMForceField" << sendl;
    this->Inherited::init();

    _topology = this->getContext()->getMeshTopology();
	this->getContext()->get(highOrderTetraGeo);

    if ((_topology->getNbTetrahedra()==0) || (!highOrderTetraGeo))
    {
        serr << "ERROR(HighOrderTetrahedralCorotationalFEMForceField): object must have a Tetrahedral Set Topology and a HighOrderTetrahedronSetGeometryAlgorithms component "<<sendl;
        return;
    }
    updateLameCoefficients();

	if (d_anisotropy.getValue() == "isotropic")
        elasticitySymmetry= ISOTROPIC;
    else if (d_anisotropy.getValue() == "cubic") 
        elasticitySymmetry= CUBIC;
    else if (d_anisotropy.getValue() == "transverseIsotropic") 
        elasticitySymmetry= TRANSVERSE_ISOTROPIC;
    else if (d_anisotropy.getValue() == "orthotropic") 
        elasticitySymmetry= ORTHOTROPIC;
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

    std::set<typename topology::NumericalIntegrationDescriptor<Real, 3>::QuadratureMethod> qmSet = highOrderTetraGeo->getTetrahedronNumericalIntegrationDescriptor().getQuadratureMethods();
    if (qmSet.count(numericalIntegrationMethod.getValue()) == 0) {
        serr << "cannot recognize numerical integration method  " << numericalIntegrationMethod.getValue() << qmSet<< sendl;
    }


    helper::vector<TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());
    tetrahedronInf.resize(_topology->getNbTetrahedra());

    if (_initialPoints.size() == 0)
    {
        // get restPosition
        const VecCoord& p = this->mstate->read(core::ConstVecCoordId::restPosition())->getValue();
        _initialPoints=p;
	}

	helper::system::thread::ctime_t startAssembly=helper::system::thread::CTime::getTime();
	totalUpdateMat=0;
	totalComputeLocalStiffness=0;

    size_t i;


	// precompute the coefficients for handling affine elements
	topology::TetrahedronIndexVector tbi1,tbi2;
	affineStiffnessCoefficientArray.clear();
	
	topology::HighOrderDegreeType degree=highOrderTetraGeo->getTopologyContainer()->getDegree();
	size_t nbControlPoints=(degree+1)*(degree+2)*(degree+3)/6;
	sofa::helper::vector<topology::TetrahedronIndexVector> tbiArray;
	tbiArray=highOrderTetraGeo->getTopologyContainer()->getTetrahedronIndexArray();

	if (degree==1) 
		integrationMethod= AFFINE_ELEMENT_INTEGRATION;

	computeElasticityTensor();
	if (degree>1) {
		if ((integrationMethod== AFFINE_ELEMENT_INTEGRATION) || (d_forceAffineAssemblyForAffineElements.getValue()))
		{
			if (degree<6) {
				affineStiffnessCoefficientPreStoredArray.allocate(degree);
			}
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
						for(m=0; m<4; ++m){
							if (m!=l) {
								coeffMatrix[l][m]-=0.5*(coeffMatrix[l][l]+coeffMatrix[m][m]);
							}
						}
					}
					if (degree==2) {
						for(l=0; l<6; ++l){
							m=edgesInTetrahedronArray[l][0];
							n=edgesInTetrahedronArray[l][1];

							(*(affineStiffnessCoefficientPreStoredArray.weightArrayQuadratic[0]))[rank][l]=coeffMatrix[m][n];
							(*(affineStiffnessCoefficientPreStoredArray.weightArrayQuadratic[1]))[rank][l]=coeffMatrix[n][m];
						}
					} else if (degree==3) {
						for(l=0; l<6; ++l){
							m=edgesInTetrahedronArray[l][0];
							n=edgesInTetrahedronArray[l][1];

							(*(affineStiffnessCoefficientPreStoredArray.weightArrayCubic[0]))[rank][l]=coeffMatrix[m][n];
							(*(affineStiffnessCoefficientPreStoredArray.weightArrayCubic[1]))[rank][l]=coeffMatrix[n][m];
						}
					} else if (degree==4) {
						for(l=0; l<6; ++l){
							m=edgesInTetrahedronArray[l][0];
							n=edgesInTetrahedronArray[l][1];

							(*(affineStiffnessCoefficientPreStoredArray.weightArrayQuartic[0]))[rank][l]=coeffMatrix[m][n];
							(*(affineStiffnessCoefficientPreStoredArray.weightArrayQuartic[1]))[rank][l]=coeffMatrix[n][m];
						}
					} else if (degree==5) {
						for(l=0; l<6; ++l){
							m=edgesInTetrahedronArray[l][0];
							n=edgesInTetrahedronArray[l][1];

							(*(affineStiffnessCoefficientPreStoredArray.weightArrayQuintic[0]))[rank][l]=coeffMatrix[m][n];
							(*(affineStiffnessCoefficientPreStoredArray.weightArrayQuintic[1]))[rank][l]=coeffMatrix[n][m];
						}

					} else {
						Vec6 coeffVec1,coeffVec2;
						for(l=0; l<6; ++l){
							m=edgesInTetrahedronArray[l][0];
							n=edgesInTetrahedronArray[l][1];
							coeffVec1[l]=coeffMatrix[m][n];
							coeffVec2[l]=coeffMatrix[n][m];
						}

						affineStiffnessCoefficientArray.push_back(coeffVec1);
						affineStiffnessCoefficientArray.push_back(coeffVec2);
					}
				}
			}
		}
		if ( (integrationMethod== NUMERICAL_INTEGRATION)  || (integrationMethod== NUMERICAL_INTEGRATION_2)  )
		{
			numericalIntegrationStiffnessDataArray.clear();
			/// get value of integration points0
			topology::NumericalIntegrationDescriptor<Real,4> &nid=highOrderTetraGeo->getTetrahedronNumericalIntegrationDescriptor();
			typename topology::NumericalIntegrationDescriptor<Real,4>::QuadraturePointArray qpa=nid.getQuadratureMethod((typename topology::NumericalIntegrationDescriptor<Real,4>::QuadratureMethod)numericalIntegrationMethod.getValue(),
				numericalIntegrationOrder.getValue());
			size_t i,j,k,l,m,n;
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
				std::vector<Vec4> shapeFunctionDerivativeArray;
				for(j=0;j<tbiArray.size();++j) {
					Vec4 deriv=highOrderTetraGeo->computeShapeFunctionDerivatives(tbiArray[j],qp.first);
					shapeFunctionDerivativeArray.push_back(deriv);
					Deriv der(deriv[0]-deriv[3],deriv[1]-deriv[3],deriv[2]-deriv[3]);
					nimd.coefficientArray.push_back(der);
				}
				size_t rank;
				for(rank=0,j=0;j<tbiArray.size();++j) {
					for(k=j+1;k<tbiArray.size();++k,++rank) {
						coeffMatrix=dyad(shapeFunctionDerivativeArray[j],shapeFunctionDerivativeArray[k])*6*weight;
						for(l=0; l<4; ++l){
							for(m=0; m<4; ++m){
								if (l!=m) {
									coeffMatrix[l][m]-=0.5*(coeffMatrix[l][l]+coeffMatrix[m][m]);
								}
							}
						}
						if (integrationMethod== NUMERICAL_INTEGRATION) { 
							if ((degree>5) || (d_oneRotationPerIntegrationPoint.getValue())) {
								Vec6 coeffVec1,coeffVec2;
								for(l=0; l<6; ++l){
									m=edgesInTetrahedronArray[l][0];
									n=edgesInTetrahedronArray[l][1];
									coeffVec1[l]=coeffMatrix[m][n];
									coeffVec2[l]=coeffMatrix[n][m];
								}
								nimd.weightArray.push_back(coeffVec1);
								nimd.weightArray.push_back(coeffVec2);
							} else {
								if (degree==2) {

									for(l=0; l<6; ++l){
										m=edgesInTetrahedronArray[l][0];
										n=edgesInTetrahedronArray[l][1];

										(*(nimd.arrayPointer.weightArrayQuadratic[0]))[rank][l]=coeffMatrix[m][n];
										(*(nimd.arrayPointer.weightArrayQuadratic[1]))[rank][l]=coeffMatrix[n][m];
									}
								} else 	if (degree==3) {

									for(l=0; l<6; ++l) {
										m=edgesInTetrahedronArray[l][0];
										n=edgesInTetrahedronArray[l][1];

										(*(nimd.arrayPointer.weightArrayCubic[0]))[rank][l]=coeffMatrix[m][n];
										(*(nimd.arrayPointer.weightArrayCubic[1]))[rank][l]=coeffMatrix[n][m];
									}
								} else 	if (degree==4) {

									for(l=0; l<6; ++l) {
										m=edgesInTetrahedronArray[l][0];
										n=edgesInTetrahedronArray[l][1];

										(*(nimd.arrayPointer.weightArrayQuartic[0]))[rank][l]=coeffMatrix[m][n];
										(*(nimd.arrayPointer.weightArrayQuartic[1]))[rank][l]=coeffMatrix[n][m];
									}
								} else 	if (degree==5) {
									for(l=0; l<6; ++l) {
										m=edgesInTetrahedronArray[l][0];
										n=edgesInTetrahedronArray[l][1];

										(*(nimd.arrayPointer.weightArrayQuintic[0]))[rank][l]=coeffMatrix[m][n];
										(*(nimd.arrayPointer.weightArrayQuintic[1]))[rank][l]=coeffMatrix[n][m];
									}
								}
							}	
						} else 	if (integrationMethod== NUMERICAL_INTEGRATION_2) { 
							nimd.weightArray4x4.push_back(coeffMatrix);
						}
					}
				}

				numericalIntegrationStiffnessDataArray.push_back(nimd);
			}
		}
        if (integrationMethod == BEZIER_NUMERICAL_INTEGRATION)
        {
            /// first fill the first vector independent from the integration points
            BezierTetrahedronSetGeometryAlgorithms<DataTypes> *bezierTetraGeo = dynamic_cast<BezierTetrahedronSetGeometryAlgorithms<DataTypes> *> (highOrderTetraGeo);
            if (bezierTetraGeo == NULL) {
                serr << "Could not find any BezierTetrahedronSetGeometryAlgorithms while using BEZIER_NUMERICAL_INTEGRATION" << sendl;
                return;
            }
            bezierCoefficientArray.clear();
            /// store the coefficient that are independent from the integration point
            topology::TetrahedronIndexVector tbi1Copy, tbi2Copy;
            size_t rank, r;
            size_t i, j, k, l, m;
            for (rank = 0, j = 0; j < tbiArray.size(); ++j) {
                for (k = j + 1; k < tbiArray.size(); ++k, ++rank) {
                    Vec16 weight;
                    tbi1 = tbiArray[j];
                    tbi2 = tbiArray[k];
                    for (r = 0, l = 0; l < 4; ++l) {
                        for (m = 0; m < 4; ++m, ++r) {
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
            for (j = 0; j < tbiArray2.size(); ++j) {
                tivMap.insert(std::make_pair(tbiArray2[j], j));
            }
            /// now fills the bezierMappingArray
            for (rank = 0, j = 0; j < tbiArray.size(); ++j) {
                for (k = j + 1; k < tbiArray.size(); ++k, ++rank) {
                    tbi1 = tbiArray[j] + tbiArray[k];
                    Vec16Int mapping;
                    for (r = 0, l = 0; l<4; ++l) {
                        for (m = 0; m<4; ++m, ++r) {

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
            topology::NumericalIntegrationDescriptor<Real, 4> &nid = highOrderTetraGeo->getTetrahedronNumericalIntegrationDescriptor();
            typename topology::NumericalIntegrationDescriptor<Real, 4>::QuadraturePointArray qpa = nid.getQuadratureMethod((typename topology::NumericalIntegrationDescriptor<Real, 4>::QuadratureMethod)numericalIntegrationMethod.getValue(),
                numericalIntegrationOrder.getValue());

            sofa::defaulttype::Vec<4, Real> bc;
            Real weight, fac;
            
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
                    nimd.weightBezierArray[j] = 6 * degree*degree*fac*bezierTetraGeo->computeShapeFunctionOfGivenDegree(tbiArray2[j], qp.first, 2 * degree - 2);
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
					Deriv der(deriv[0]-deriv[3],deriv[1]-deriv[3],deriv[2]-deriv[3]);
					nimd.coefficientArray.push_back(der);
				}
				numericalIntegrationStiffnessDataArray.push_back(nimd);
			}



		}
		helper::system::thread::ctime_t startComputeLocalStiffness=helper::system::thread::CTime::getTime();
		/// initialize the data structure associated with each tetrahedron
	}
    for (i=0; i<_topology->getNbTetrahedra(); ++i)
    {
        tetrahedronHandler->applyCreateFunction(i,tetrahedronInf[i],_topology->getTetrahedron(i),
                (const helper::vector< unsigned int > )0,
                (const helper::vector< double >)0);
    }

    //msg_info() << getStiffnessArray(0,&tetrahedronInf[0]);

	helper::system::thread::ctime_t endComputeLocalStiffness=helper::system::thread::CTime::getTime();
	if (this->f_printLog.getValue()) {
		helper::system::thread::ctime_t endAssembly=helper::system::thread::CTime::getTime();
		std::cerr<< "Assembly time="<< ((endAssembly-startAssembly)/(Real)helper::system::thread::CTime::getRefTicksPerSec()) << std::endl;
		std::cerr<<" total update mat="<<((totalUpdateMat)/(Real)helper::system::thread::CTime::getRefTicksPerSec()) << std::endl;
		
		std::cerr<<" total compute local stiffness ="<<((totalComputeLocalStiffness)/(Real)helper::system::thread::CTime::getRefTicksPerSec()) << std::endl;
		
		d_assemblyTime.setValue(((endAssembly-startAssembly)/(Real)helper::system::thread::CTime::getRefTicksPerSec()));

		//startAssembly=helper::system::thread::CTime::getTime();
		//helper::system::thread::CTime::sleep(1.0);
		//endAssembly=helper::system::thread::CTime::getTime();
		//std::cerr<< "test time="<< ((endAssembly-startAssembly)/(Real)helper::system::thread::CTime::getRefTicksPerSec()) << std::endl;
	}
    /// set the call back function upon creation of a tetrahedron
    tetrahedronInfo.createTopologicalEngine(_topology,tetrahedronHandler);
    tetrahedronInfo.registerTopologicalData();
    tetrahedronInfo.endEdit();

	updateTopologyInfo=true;

}


template <class DataTypes>
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::updateTopologyInformation()
{
    int i;
    unsigned int j;

    int nbTetrahedra=_topology->getNbTetrahedra();

    TetrahedronRestInformation *tetinfo;

    helper::vector<typename HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::TetrahedronRestInformation>& tetrahedronInf = *(tetrahedronInfo.beginEdit());
   


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
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::computeElasticityTensor() 											
{
	const helper::vector<Real> & anisotropyParameter=d_anisotropyParameter.getValue();
	const helper::vector<Coord> & anisotropyDirection=d_anisotropyDirection.getValue();

	if (elasticitySymmetry==ISOTROPIC) {
        // elasticity tensor in isotropic case
        // lambda = E*v/(1-2*v)*(1+v)
        // mu = E/(2*(1+v))
        // with v == poissonRatio & E youngModulus
        Real lambda=getLambda();
        Real mu=getMu();
        //
        //                        +---+-----+-----+-----+----------+----------+----------+
        //                        | / | 0   | 1   | 2   | 3        | 4        | 5        |
        //                        +===+=====+=====+=====+==========+==========+==========+
        //                        | 0 | 1-v | v   | v   | 0        | 0        | 0        |
        //                        +---+-----+-----+-----+----------+----------+----------+
        //                        | 1 | v   | 1-v | v   | 0        | 0        | 0        |
        //                        +---+-----+-----+-----+----------+----------+----------+
        //  E/(1+v)(1-2v) *       | 2 | v   | v   | 1-v | 0        | 0        | 0        |
        //                        +---+-----+-----+-----+----------+----------+----------+
        //                        | 3 | 0   | 0   | 0   | (1-2v)/2 | 0        | 0        |
        //                        +---+-----+-----+-----+----------+----------+----------+
        //                        | 4 | 0   | 0   | 0   | 0        | (1-2v)/2 | 0        |
        //                        +---+-----+-----+-----+----------+----------+----------+
        //                        | 5 | 0   | 0   | 0   | 0        | 0        | (1-2v)/2 |
        //                        +---+-----+-----+-----+----------+----------+----------+
        //
        //     ==
        //                        +---+-------------+-------------+-------------+----+----+----+
        //                        | / | 0           | 1           | 2           | 3  | 4  | 5  |
        //                        +===+=============+=============+=============+====+====+====+
        //                        | 0 | 2*mu+lambda | lambda      | lambda      | 0  | 0  | 0  |
        //                        +---+-------------+-------------+-------------+----+----+----+
        //                        | 1 | lambda      | 2*mu+lambda | lambda      | 0  | 0  | 0  |
        //                        +---+-------------+-------------+-------------+----+----+----+
        //                        | 2 | lambda      | lambda      | 2*mu+lambda | 0  | 0  | 0  |
        //                        +---+-------------+-------------+-------------+----+----+----+
        //                        | 3 | 0           | 0           | 0           | mu | 0  | 0  |
        //                        +---+-------------+-------------+-------------+----+----+----+
        //                        | 4 | 0           | 0           | 0           | 0  | mu | 0  |
        //                        +---+-------------+-------------+-------------+----+----+----+
        //                        | 5 | 0           | 0           | 0           | 0  | 0  | mu |
        //                        +---+-------------+-------------+-------------+----+----+----+

        elasticityTensor(0,0)=2*mu+lambda;
        elasticityTensor(1,1)=2*mu+lambda;
        elasticityTensor(2,2)=2*mu+lambda;
        elasticityTensor(3,3)=mu;
        elasticityTensor(4,4)=mu;
        elasticityTensor(5,5)=mu;

		elasticityTensor(0,1)=lambda;elasticityTensor(0,2)=lambda;elasticityTensor(1,2)=lambda;
		elasticityTensor(1,0)=lambda;elasticityTensor(2,0)=lambda;elasticityTensor(2,1)=lambda;

        //msg_info() << "C11 "<< elasticityTensor(0,0);
        //msg_info() << "C12 "<< elasticityTensor(1,0);
        //msg_info() << "C44 "<< elasticityTensor(4,4);

    }
    else
    {
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

        // Here we will use the work of Sandrine Germian see thesis (2015): https://opus4.kobv.de/opus4-fau/frontdoor/index/index/docId/3490
        // Chap 4 Spectral decomposition and the Kelvin modes :
        //
        // - For CUBIC                  -------> 4.4.1 The cubic crystal system Materials (p38)
        //
        // - For TRANSVERSE_ISOTROPIC   -------> 4.4.6 The tetragonal crystal system Materials (p48)
        //
        // As well as the equations from http://solidmechanics.org/text/Chapter3_2/Chapter3_2.htm


        if (elasticitySymmetry==CUBIC) {

            // get the different constants : young modulus, Poisson ratio and anisotropy ratio.
            Real youngModulus=d_youngModulus.getValue();
            Real poissonRatio=d_poissonRatio.getValue();
            Real anisotropyRatio=anisotropyParameter[0];

            // Same as the isotropic but with a coefficient of anysotropy (we will call A)
            // if A == 1    --> the material is isotropic
            // defined by the expression :
            //
            // mu = A * E/(2*(1+v))
            //
            //
            //                        +---+-----+-----+-----+----------+----------+-----------+
            //                        | / | 0   | 1   | 2   | 3        | 4        | 5         |
            //                        +===+=====+=====+=====+==========+==========+===========+
            //                        | 0 | 1-v | v   | v   | 0        | 0        | 0         |
            //                        +---+-----+-----+-----+----------+----------+-----------+
            //                        | 1 | v   | 1-v | v   | 0        | 0        | 0         |
            //                        +---+-----+-----+-----+----------+----------+-----------+
            //  E/(1+v)(1-2v) *       | 2 | v   | v   | 1-v | 0        | 0        | 0         |
            //                        +---+-----+-----+-----+----------+----------+-----------+
            //                        | 3 | 0   | 0   | 0   | A(1-2v)/2 | 0        | 0        |
            //                        +---+-----+-----+-----+----------+----------+-----------+
            //                        | 4 | 0   | 0   | 0   | 0        | A(1-2v)/2 | 0        |
            //                        +---+-----+-----+-----+----------+----------+-----------+
            //                        | 5 | 0   | 0   | 0   | 0        | 0        | A(1-2v)/2 |
            //                        +---+-----+-----+-----+----------+----------+-----------+

            Real c11 = youngModulus * (1 - poissonRatio) / (1 - poissonRatio - 2 * poissonRatio * poissonRatio);
            Real c12 = youngModulus * poissonRatio / (1 - poissonRatio - 2 * poissonRatio * poissonRatio);
            Real c44 = anisotropyRatio * (c11 - c12) / 2;

            // The number of modes is equal to three, so that the tensor has three eigenvalues:
            Real eigen1 = c11 + 2 * c12;    // (4.39) dim 1
            Real eigen2 = c11 - c12;        // (4.40) dim 2
            Real eigen3 = 2 * c44;          // (4.41) dim 3

            //msg_info() << "anisotropyRatio "<<anisotropyRatio ;
            //msg_info() << "C11 "<< c11;
            //msg_info() << "C12 "<< c12;
            //msg_info() << "C44 (mu) "<< c44;

            // ----------------------------------------------------------------------
            // BUILD THE ORTHOGONAL MATRICES
            // ----------------------------------------------------------------------

            // The spectral decomposition of the cubic tensor consists of a dilatation,
            // a two-dimensional and a three-dimensional eigenspace
            //  (here 'x' is the dyadic product)

            // The first projection tensor associated with the first eigenvalue depends on the dilatation mode as for the isotropic tensor
            // P1 = Nd x Nd

            Mat3x3 Nd;
            Nd.identity();
            Nd/= sqrt(3.0f);

            // The second projection tensor associated with the second eigenvalue is the sum of the tensor product of the isochoric extension and pure shear modes with themselves
            // P2 = Nei x Nei + Npi x Npi (with i == 1/2/3)

            Mat3x3 Ne=(2*dyad(n,n)-dyad(v1,v1)-dyad(v2,v2))/sqrt(6.0f);
            Mat3x3 Np=(dyad(v1,v1)-dyad(v2,v2))/sqrt(2.0f);

            // The third projection tensor associated with the third eigenvalue is the sum of the tensor product of the three simple shear modes with themselves
            // P3 = Ns1 x Ns1 + Ns2 x Ns2 + Ns3 x Ns3

            Mat3x3 Ns1=(dyad(v2,n)+dyad(n,v2))/sqrt(2.0f);
            Mat3x3 Ns2=(dyad(v1,n)+dyad(n,v1))/sqrt(2.0f);
            Mat3x3 Ns3=(dyad(v1,v2)+dyad(v2,v1))/sqrt(2.0f);


            // push all symmetric matrices and the eigenvalues
            anisotropyMatrixArray.push_back(Nd);
            anisotropyScalarArray.push_back(eigen1);
            anisotropyMatrixArray.push_back(Ne);
            anisotropyScalarArray.push_back(eigen2);
            anisotropyMatrixArray.push_back(Np);
            anisotropyScalarArray.push_back(eigen2);
            anisotropyMatrixArray.push_back(Ns1);
            anisotropyScalarArray.push_back(eigen3);
            anisotropyMatrixArray.push_back(Ns2);
            anisotropyScalarArray.push_back(eigen3);
            anisotropyMatrixArray.push_back(Ns3);
            anisotropyScalarArray.push_back(eigen3);
        }
        else if (elasticitySymmetry==TRANSVERSE_ISOTROPIC) {

            // get the constants from the young modulus, Poisson ratio and anisotropy ratio.
            Real youngModulusTransverse=d_youngModulus.getValue();
            Real poissonRatioTransverse=d_poissonRatio.getValue();
            Real youngModulusLongitudinal=anisotropyParameter[0];
            Real poissonRatioTransverseLongitudinal=anisotropyParameter[1];
            Real shearModulusTransverse=anisotropyParameter[2];

            Real poissonRatioLongitudinalTransverse=poissonRatioTransverseLongitudinal*youngModulusLongitudinal/youngModulusTransverse;
            Real gamma=1/(1-poissonRatioTransverse*poissonRatioTransverse-2*poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal-2*poissonRatioTransverse*poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal);

            Real c11=youngModulusTransverse*gamma*(1-poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal);
            Real c33=youngModulusLongitudinal*gamma*(1-poissonRatioTransverse*poissonRatioTransverse);
            Real c12=youngModulusTransverse*gamma*(poissonRatioTransverse+poissonRatioLongitudinalTransverse*poissonRatioTransverseLongitudinal);
            Real c13=youngModulusTransverse*gamma*(poissonRatioLongitudinalTransverse+poissonRatioTransverse*poissonRatioTransverseLongitudinal);
            Real c66=youngModulusTransverse/(2*(1+poissonRatioTransverse));
            Real c44=shearModulusTransverse; // = c55

            //msg_info() << "C11      "<< c11;
            //msg_info() << "C12      "<< c12;
            //msg_info() << "C33      "<< c33;
            //msg_info() << "C44/C55  "<< c44;
            //msg_info() << "C66      "<< c66;

            // (4.119)
            Real talpha=sqrt(2.0f)*(c11+c12-c33)/(4*c13);
            Real alpha=atan(talpha);
            Real salpha=sin(alpha);
            Real calpha=cos(alpha);
            Real secalpha=1/calpha;

            // (4.118)
            Real eigen1 = c33 + M_SQRT2 * c13 * (talpha + secalpha);
            Real eigen2 = c11 - c12;
            Real eigen3 = c33 + M_SQRT2 * c13 * (talpha - secalpha);
            Real eigen4 = 2 * shearModulusTransverse;

            // ----------------------------------------------------------------------
            // BUILD THE ORTHOGONAL MATRICES
            // ----------------------------------------------------------------------

            // No dilatation/extension Kelvin modes ? Nd/Ne
            // yes see p48 :
            // because we represent the tetragonal tensor in the e3-direction
            // The five projection tensors are expressed as functions of the typical Kelvin modes and two dilatation

            // (4.120)
            //P1 = Nh1 x Nh1
            //P2 = Np3 x Np3 + Ns3 x Ns3
            //P3 = Nh2 x Nh2
            //P4 = Ns1 x Ns1 + Ns2 x Ns2

            // ---------------------------------------------------------
            // dilatation modes
            // defined in Equation 4.121 & 4.122 (p51)
            Real val1_Nh1 = 0.5*(1+salpha)+sqrt(2.0)*calpha/4.0f;
            Real val2_Nh1 = 0.5*(1-salpha)+sqrt(2.0)*calpha/2.0f;
            Real val1_Nh2 = 0.5*(1-salpha)-sqrt(2.0)*calpha/4.0f;
            Real val2_Nh2 = 0.5*(1+salpha)-sqrt(2.0)*calpha/2.0f;

            Mat3x3 Nh1 = val1_Nh1 * ( dyad(v1,v1) + dyad(v2,v2) ) + val2_Nh1 * dyad(n,n);
            Mat3x3 Nh2 = val1_Nh2 * ( dyad(v1,v1) + dyad(v2,v2) ) + val2_Nh2 * dyad(n,n);

            // normalization Nh1/Nh2
            Nh1/=sqrt(2*val1_Nh1*val1_Nh1+val2_Nh1*val2_Nh1);
            Nh2/=sqrt(2*val1_Nh2*val1_Nh2+val2_Nh2*val2_Nh2);

            // ---------------------------------------------------------
            // isochoric pure shear modes along e3
            // defined in Equation 4.18 (p35)
            Mat3x3 Np=(dyad(v1,v1)-dyad(v2,v2))/sqrt(2.0f);

            // ---------------------------------------------------------
            // The three isochoric simple shear modes along e1, e2, and e3
            // defined in Equation 4.19, Equation 4.20, and Equation 4.21, respectively (p35)
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
}


template<class DataTypes>
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::computeTetrahedronStiffnessEdgeMatrix(const Coord point[4], Mat6x9 edgeStiffnessVectorized[2])
{
	helper::system::thread::ctime_t startComputeStiffness=helper::system::thread::CTime::getTime();

	Coord shapeVector[4];
	Mat3x3 edgeStiffness[6];
	/// compute 6 times the rest volume
	Real volume=dot(cross(point[1]-point[0],point[2]-point[0]),point[0]-point[3]);
	/// store the rest volume
    // my_tinfo.restVolume=volume/6;

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
                edgeStiffness[j]+=(anisotropyScalarArray[i]*anisotropyMatrixArray[i]*tmp*anisotropyMatrixArray[i])*fabs(volume)/6;
			}
		}
	}
	size_t p;
	for(j=0; j<6; ++j)
	{
		k=edgesInTetrahedronArray[j][0];
		l=edgesInTetrahedronArray[j][1];
		for(p=0,m=0; m<3; ++m)
		{
			for(n=0; n<3; ++n,++p)
			{
				edgeStiffnessVectorized[0][j][p]=edgeStiffness[j][m][n];
				edgeStiffnessVectorized[1][j][p]=edgeStiffness[j][n][m];
			}
		}
	}
	if (this->f_printLog.getValue()) {
		helper::system::thread::ctime_t endComputeStiffness=helper::system::thread::CTime::getTime();
		totalComputeLocalStiffness+=endComputeStiffness-startComputeStiffness;
	}
}


template<class DataTypes>
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::computeTetrahedronStiffnessEdgeMatrix(const Coord point[4], Mat3x3 edgeStiffness[6])
{
	helper::system::thread::ctime_t startComputeStiffness=helper::system::thread::CTime::getTime();

	Coord shapeVector[4];
	/// compute 6 times the rest volume
	Real volume=dot(cross(point[1]-point[0],point[2]-point[0]),point[0]-point[3]);
	/// store the rest volume
    // my_tinfo.restVolume=volume/6;

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
                edgeStiffness[j]+=(anisotropyScalarArray[i]*anisotropyMatrixArray[j]*tmp*anisotropyMatrixArray[j])*fabs(volume)/6;
			}
		}
	}
	if (this->f_printLog.getValue()) {
		helper::system::thread::ctime_t endComputeStiffness=helper::system::thread::CTime::getTime();
		totalComputeLocalStiffness+=endComputeStiffness-startComputeStiffness;
	}

}


template<class DataTypes>
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::computeQRRotation( Mat3x3 &r, const Coord *dp)
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
const helper::vector<typename HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::Mat3x3> &
    HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::getStiffnessArray(
    const size_t i,
    typename HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::TetrahedronRestInformation *restTetra)
{
    return(restTetra->stiffnessVector);
}

template <class DataTypes>
const  helper::vector<typename HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::Mat3x3> &
    HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::getRotatedStiffnessArray(
    const size_t i,
    const typename HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::TetrahedronRestInformation *restTetra)
{

    if (decompositionMethod==LINEAR_ELASTIC)
    {
        return(restTetra->stiffnessVector);
    } else
        return(restTetra->rotatedStiffnessVector);
}

template <class DataTypes>
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* mparams /* PARAMS FIRST */, 
																	 DataVecDeriv &  dataF, const DataVecCoord &  dataX , const DataVecDeriv & /*dataV*/ )
{

    VecDeriv& f        = *(dataF.beginEdit());
    const VecCoord& x  =   dataX.getValue()  ;
    const VecCoord& x0= this->mstate->read(core::ConstVecCoordId::restPosition())->getValue();
	sofa::helper::vector<Coord> dp,force;
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
	
	dp.resize(nbControlPoints);
	force.resize(nbControlPoints);
	
	Coord dpos,sv;
		
	for(i=0; i<nbTetrahedra; i++ )
	{
		tetinfo=&tetrahedronInf[i];
		const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=
			highOrderTetraGeo->getTopologyContainer()->getGlobalIndexArrayOfControlPoints(i);
		
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
				highOrderTetraGeo->computeNodalValueDerivatives(i,numericalIntegrationStiffnessDataArray[l].integrationPoint, x,point);

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
				const std::vector<Vec6> & weightArray=numericalIntegrationStiffnessDataArray[l].weightArray;
				const Mat6x9  &edgeStiff1 =tetinfo->integrationPointsStiffnessVector[2*l];
				const Mat6x9  &edgeStiff2 =tetinfo->integrationPointsStiffnessVector[2*l+1];

				for (rank=0,j=0; j<nbControlPoints; ++j) {
					v0 = indexArray[j];
					for ( k=j+1; k<nbControlPoints; ++k,++rank) {

						const Vec6  & coeffVec1=weightArray[2*rank];
						const Vec6  & coeffVec2=weightArray[2*rank+1];
						Vec9 res=edgeStiff1.multTranspose(coeffVec1)+edgeStiff2.multTranspose(coeffVec2);
						Mat3x3 stiffness= Mat3x3((const Real *) &res[0]);
						//			Vec9 res=edgeStiffness[0]*coeffVec1+edgeStiffness[1]*coeffVec2;


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
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::addDForce(const sofa::core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv&   datadF , const DataVecDeriv&   datadX )
{
    VecDeriv& df       = *(datadF.beginEdit());
    const VecCoord& dx =   datadX.getValue()  ;
    Real kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());
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
		highOrderTetraGeo->getTopologyContainer()->getGlobalIndexArrayOfHighOrderPointsInTetrahedron(i,indexArray) ;
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
SReal HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::getPotentialEnergy(const core::MechanicalParams*, const DataVecCoord&) const
{
    serr << "ERROR("<<this->getClassName()<<"): getPotentialEnergy( const MechanicalParams*, const DataVecCoord& ) not implemented." << sendl;
    return 0.0;
}


template<class DataTypes>
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::updateLameCoefficients()
{
    lambda= d_youngModulus.getValue()*d_poissonRatio.getValue()/((1-2*d_poissonRatio.getValue())*(1+d_poissonRatio.getValue())); // E*v/(1-2*v)*(1+v)
    mu = d_youngModulus.getValue()/(2*(1+d_poissonRatio.getValue())); // E/(2*(1+v))

}


template<class DataTypes>
void HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
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

#endif // SOFA_COMPONENT_FORCEFIELD_HighOrderTetrahedralCorotationalFEMForceField_INL
