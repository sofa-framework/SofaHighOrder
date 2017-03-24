#ifndef SOFA_COMPONENT_FORCEFIELD_ADAPTIVEBEZIERTETRAHEDRALCOROTATIONALFEMFORCEFIELD_INL
#define SOFA_COMPONENT_FORCEFIELD_ADAPTIVEBEZIERTETRAHEDRALCOROTATIONALFEMFORCEFIELD_INL

#include "AdaptiveBezierTetrahedralCorotationalFEMForceField.h"
#include "AdaptiveBezierTetrahedronSetTopologyContainer.h"
#include <sofa/core/visual/VisualParams.h>
#include <BezierTetrahedronSetGeometryAlgorithms.h>
#include "AdaptiveBezierTetrahedronSetTopologyModifier.h"

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
typedef sofa::defaulttype::Vec<4,unsigned int> Vec4ui;

const unsigned int edgesInTetrahedronArray[6][2] = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}};
     
template <class DataTypes> AdaptiveBezierTetrahedralCorotationalFEMForceField<DataTypes>::AdaptiveBezierTetrahedralCorotationalFEMForceField() :
  f_indexType(initData(&f_indexType, std::string("strain"),"indexType", "the type of Index used to measure the deformation of a Bezier tetrahedron"))
, f_deformationIndexArray(initData(&f_deformationIndexArray, "deformationIndex", "an array of scalar describing the amount of deformation of a tetrahedron"))

{
   
}

template <class DataTypes> AdaptiveBezierTetrahedralCorotationalFEMForceField<DataTypes>::~AdaptiveBezierTetrahedralCorotationalFEMForceField()
{
   
}

template <class DataTypes> void AdaptiveBezierTetrahedralCorotationalFEMForceField<DataTypes>::init()
{
    //	serr << "initializing HighOrderTetrahedralCorotationalFEMForceField" << sendl;
    this->HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::init();

	this->getContext()->get(adaptiveContainer);

    if ((this->_topology->getNbTetrahedra()==0) || (!adaptiveContainer))
    {
        serr << "ERROR(AdaptiveBezierTetrahedralCorotationalFEMForceField): object must have a AdaptiveBezierTetrahedronSetTopologyContainer  component "<<sendl;
        return;
    }
	// initialize reducedStiffnessVector as a copy of stiffnessVector
    helper::vector<typename HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::TetrahedronRestInformation>& tetrahedronInf = *(this->tetrahedronInfo.beginEdit());
     for (size_t i=0; i<this->_topology->getNbTetrahedra(); ++i)
    {
		tetrahedronInf[i].reducedStiffnessVector.resize(tetrahedronInf[i].stiffnessVector.size());
		std::copy(tetrahedronInf[i].stiffnessVector.begin(),tetrahedronInf[i].stiffnessVector.end(),
			tetrahedronInf[i].reducedStiffnessVector.begin());
	 }
}

/*
template <class DataTypes>
void AdaptiveBezierTetrahedralCorotationalFEMForceField<DataTypes>::updateTopologyInformation()
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
*/
template<class DataTypes>
void AdaptiveBezierTetrahedralCorotationalFEMForceField<DataTypes>::handleTopologyChange(core::topology::Topology *topo)
{
	core::topology::BaseMeshTopology*	_topology = dynamic_cast<core::topology::BaseMeshTopology*> (topo);

	if(_topology != NULL)
	{
		std::list<const core::topology::TopologyChange *>::const_iterator itBegin=_topology->beginChange();
		std::list<const core::topology::TopologyChange *>::const_iterator itEnd=_topology->endChange();

		while( itBegin != itEnd )
		{
			if ((*itBegin)->getChangeType()== core::topology::TOPOLOGYCHANGE_LASTID) {
				const HighOrderTetrahedronDegreeChanged *btdc=dynamic_cast<const HighOrderTetrahedronDegreeChanged *>(*itBegin);
					if (btdc) {
                        typename HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::TetrahedronRestInformation *tetinfo;
						size_t nbDofs;
                        HighOrderDegreeType degree=this->highOrderTetraGeo->getTopologyContainer()->getDegree();
						size_t nbControlPoints=(degree+1)*(degree+2)*(degree+3)/6;
                        helper::vector<typename HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::TetrahedronRestInformation>& tetrahedronInf = *(this->tetrahedronInfo.beginEdit());
						const sofa::helper::vector<unsigned int> &tetraModified=btdc->getTetrahedronArray();
						size_t j,k,l,m,n,p,q,index;
						AdaptiveHighOrderTetrahedronSetTopologyContainer::WeightedDOFArray wda[2];
						AdaptiveHighOrderTetrahedronSetTopologyContainer::WeightedDOF wd[2];
						
						for (size_t i=0;i<tetraModified.size();++i) {
							// updates the 
							tetinfo=&tetrahedronInf[tetraModified[i]];
							const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=
							adaptiveContainer->getTetrahedronDOFArray(tetraModified[i]);
							const AdaptiveHighOrderTetrahedronSetTopologyContainer::SeqWeightedDOFArray &swda=
								adaptiveContainer->getWeightedDOFArray(tetraModified[i]);
							nbDofs=indexArray.size();
							// resize and fill with zero entries
							tetinfo->reducedStiffnessVector.resize(nbDofs*(nbDofs-1)/2);
							std::fill(tetinfo->reducedStiffnessVector.begin(),
								tetinfo->reducedStiffnessVector.end(),Mat3x3());
							std::vector<Mat3x3> diagonalStifffnessVector; // store diagonal terms
							diagonalStifffnessVector.resize(nbControlPoints);
							std::fill(diagonalStifffnessVector.begin(),
								diagonalStifffnessVector.end(),Mat3x3());
							// handle offdiagonal stiffness terms
							for ( l=0,j=0;j<nbControlPoints;++j) {
								wda[0]=swda[j];
								for ( k=j+1;k<nbControlPoints;++k,++l) {
									wda[1]=swda[k];
									// stiffness matrix to be split among the mapped DOF
									const Mat3x3 &mat=tetinfo->stiffnessVector[l];
									// update diagonal terms as the opposite sum of diagonal terms
									diagonalStifffnessVector[k]-=mat;
									diagonalStifffnessVector[j]-=mat.transposed();
									for (m=0;m<wda[0].size();++m) {
										wd[0]=wda[0][m];
										for (n=0;n<wda[1].size();++n) {
											wd[1]=wda[1][n];
											// only consider off-diagonal term ignoring diagonal terms
											if (wd[0].first!=wd[1].first) {
											//	TetrahedronIndexVector tbi[2];
											//	tbi[0]=adaptiveContainer->getTetrahedronIndexVector(j);
											//	tbi[1]=adaptiveContainer->getTetrahedronIndexVector(k);
												// formula to get the index of the offdiagonal term in the array of stiffnesses
												if (wda[0][m].first<wda[1][n].first) {
													index=wd[0].first*nbDofs-wd[0].first*(wd[0].first+1)/2+(wd[1].first-wd[0].first-1);
													tetinfo->reducedStiffnessVector[index]+=mat*wd[0].second*wd[1].second;
												} else {
													index=wd[1].first*nbDofs-wd[1].first*(wd[1].first+1)/2+(wd[0].first-wd[1].first-1);
													tetinfo->reducedStiffnessVector[index]+=mat.transposed()*wd[0].second*wd[1].second;
												}
												/*
												std::cerr<<" input vertex "<<j<<" with tbi= "<<(Vec4ui)tbi[0]<< " and "<< k<<" with tbi= "<<(Vec4ui)tbi[1]<< " with edge stiffness rank "<<l<<" with mat= "<<mat<<std::endl;
												
												std::cerr<<" adding to reduced stiffness with index  "<<index<<" with weights "<<wd[0].second*wd[1].second<<std::endl;
												std::cerr << " local reduced vertex 0 = "<<wd[0].first << " local reduced vertex 1 = "<<wd[1].first<<std::endl;
												std::cerr << "global reduced vertex 0 = "<<indexArray[wd[0].first] << "global  reduced vertex 1 = "<<indexArray[wd[1].first]<<std::endl;
												std::cerr << " weight reduced vertex 0 = "<<wd[0].second << "  weight reduced vertex 1 = "<<wd[1].second<<std::endl; */
											}
										}
									}
								}
							}
//							for ( j=0;j<nbDofs;++j) {
//								std::cerr<<" reduced vertex index["<<j<<"]= "<<indexArray[j]<<std::endl;
//							}
							// handle diagonal terms
							for ( j=0;j<nbControlPoints;++j) {
								wda[0]=swda[j];

								// stiffness matrix to be split among the mapped DOF
								const Mat3x3 &mat=diagonalStifffnessVector[j];
								assert(fabs(mat[1][2]-mat[2][1])<1e-6);
								if (wda[0].size()>1) {
									for (m=0;m<wda[0].size();++m) {
										for (n=m+1;n<wda[0].size();++n) {
											// sort the indidces of the 2 weighted array
											if (wda[0][m].first<wda[0][n].first) {
												wd[0]=wda[0][m];
												wd[1]=wda[0][n];
											} else {
												wd[1]=wda[0][m];
												wd[0]=wda[0][n];
											}
										
											// formula to get the index of the offdiagonal term in the array of stiffnesses
											index=wd[0].first*nbDofs-wd[0].first*(wd[0].first+1)/2+(wd[1].first-wd[0].first-1);
											tetinfo->reducedStiffnessVector[index]+=mat*wd[0].second*wd[1].second;

										}
									}
								}
							}
	/*						for ( l=0,j=0;j<nbDofs;++j) {
								
								for ( k=j+1;k<nbDofs;++k,++l) {
										std::cerr<<" reduced stiffness["<<l<<"]= "<<tetinfo->reducedStiffnessVector[l]<<std::endl;
								}
							}*/
#ifndef NDEBUG
							// check the rank of each off diagnoal stiffness tensor
							for ( l=0,j=0;j<nbDofs;++j) {
								
								for ( k=j+1;k<nbDofs;++k,++l) {
	//								assert(fabs(determinant(tetinfo->reducedStiffnessVector[l]))>1e-7);
								}
							}
#endif
						}
					}
			}
			++itBegin;
		}
	}
}

template <class DataTypes>
const helper::vector<typename AdaptiveBezierTetrahedralCorotationalFEMForceField<DataTypes>::Mat3x3> &
	AdaptiveBezierTetrahedralCorotationalFEMForceField<DataTypes>::getStiffnessArray(
	const size_t i,
    typename HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::TetrahedronRestInformation *restTetra)
{

	
	//this->bezierTetraGeo->getTopologyContainer()->getGlobalIndexArrayOfBezierPointsInTetrahedron(i,indexArray) ;
	return(restTetra->reducedStiffnessVector);
}

template <class DataTypes>
const helper::vector<typename AdaptiveBezierTetrahedralCorotationalFEMForceField<DataTypes>::Mat3x3> &
	AdaptiveBezierTetrahedralCorotationalFEMForceField<DataTypes>::getRotatedStiffnessArray(
	const size_t i,
    const typename HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::TetrahedronRestInformation *restTetra)
{
    if (this->decompositionMethod==this->LINEAR_ELASTIC)
	{
		return(restTetra->reducedStiffnessVector);
	} else
		return(restTetra->rotatedStiffnessVector);
}
template <class DataTypes>
void AdaptiveBezierTetrahedralCorotationalFEMForceField<DataTypes>::update() {
}

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_HighOrderTetrahedralCorotationalFEMForceField_INL
