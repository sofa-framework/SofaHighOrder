
#ifndef SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERMESHMATRIXMASS_INL
#define SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERMESHMATRIXMASS_INL

#include "AdaptiveBezierMeshMatrixMass.h"
#include "AdaptiveBezierTetrahedronSetTopologyContainer.h"
#include "BezierTetrahedronSetGeometryAlgorithms.h"
#include "AdaptiveBezierTetrahedronSetTopologyModifier.h"
#include <SofaBaseTopology/TopologyData.inl>


namespace sofa
{

namespace component
{

namespace mass
{

template <class DataTypes, class MassType>
AdaptiveBezierHighOrderMeshMatrixMass<DataTypes, MassType>::AdaptiveBezierHighOrderMeshMatrixMass() :
 reducedTetrahedronMassInfo( initData(&reducedTetrahedronMassInfo, "reducedTetrahedronMass","values of the particles masses for all control points inside a Bezier tetrahedron") )

{
 
}


template <class DataTypes, class MassType>
AdaptiveBezierHighOrderMeshMatrixMass<DataTypes, MassType>::~AdaptiveBezierHighOrderMeshMatrixMass() 
{
 
}


template <class DataTypes, class MassType>
void AdaptiveBezierHighOrderMeshMatrixMass<DataTypes, MassType>::init()
{
   

    this->Inherited::init();
   

	this->getContext()->get(adaptiveContainer);
    if ((this->_topology->getNbTetrahedra()==0) || (!adaptiveContainer))
	{
		serr << "ERROR(AdaptiveHighOrderTetrahedralCorotationalFEMForceField): object must have a AdaptiveBezierTetrahedronSetTopologyContainer  component "<<sendl;
		return;
	}
	// initialize reducedStiffnessVector as a copy of stiffnessVector
    const MassVectorVector& tetrahedronMass = this->tetrahedronMassInfo.getValue();
	MassVectorVector& reducedTetrahedronMass = *reducedTetrahedronMassInfo.beginEdit();

	reducedTetrahedronMass.resize(tetrahedronMass.size());
    for (size_t i=0; i<this->_topology->getNbTetrahedra(); ++i)
	{
		reducedTetrahedronMass[i].resize(tetrahedronMass[i].size());
		std::copy(tetrahedronMass[i].begin(),tetrahedronMass[i].end(),
			reducedTetrahedronMass[i].begin());
	}

	reducedTetrahedronMassInfo.endEdit();
	checkConsistencyOfLumpedMass();
}

template <class DataTypes, class MassType>
void AdaptiveBezierHighOrderMeshMatrixMass<DataTypes, MassType>::checkConsistencyOfLumpedMass()
{
#ifndef NDEBUG
	// check consistency between the mass lumped array vertexMassInfo and  Bezier mass stored in reducedTetrahedronMassInfo
	size_t i,j,k,l;
	HighOrderTetrahedronSetTopologyContainer::HighOrderTetrahedronPointLocation location;
	size_t elementIndex,elementOffset,tetrahedronIndex,index;
	Real error;
    const helper::vector<MassType> &vertexMass = this->vertexMassInfo.getValue();
	const MassVectorVector& reducedTetrahedronMass = reducedTetrahedronMassInfo.getValue();
	MassType mass,totalMass;
	totalMass=(MassType)0.0f;
	for (i=0;i<adaptiveContainer->getNbPoints();++i) {
		mass=(MassType)0.0f;
		adaptiveContainer->getLocationFromGlobalIndex(i,location,elementIndex,elementOffset);
		if (location==HighOrderTetrahedronSetTopologyContainer::TETRAHEDRON) {
			// finds the mass associated with that point in the tetrahedron
			const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=
				adaptiveContainer->getTetrahedronDOFArray(elementIndex);
			for (j=0;indexArray[j]!=i;++j);
			assert(j<indexArray.size());
			// now sum the mass for each row element of the  local mass matrix
			for (k=0;k<indexArray.size();++k) {
				if (k<j) {
					index=k*indexArray.size()-k*(k-1)/2+j-k;

				} else {
					index=j*indexArray.size()-j*(j-1)/2+k-j;
				}
				mass+=reducedTetrahedronMass[elementIndex][index];
			}
		} else if (location==HighOrderTetrahedronSetTopologyContainer::TRIANGLE) {
			const TetrahedraAroundTriangle &tat=adaptiveContainer->getTetrahedraAroundTriangle(elementIndex);
			for (l=0;l<tat.size();++l) {
				tetrahedronIndex=tat[l];
				// finds the mass associated with that point in the tetrahedron
				const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=
					adaptiveContainer->getTetrahedronDOFArray(tetrahedronIndex);
				for (j=0;indexArray[j]!=i;++j);
				assert(j<indexArray.size());
				// now sum the mass for each row element of the  local mass matrix
				for (k=0;k<indexArray.size();++k) {
					if (k<j) {
						index=k*indexArray.size()-k*(k-1)/2+j-k;
					} else {
						index=j*indexArray.size()-j*(j-1)/2+k-j;
					}
					mass+=reducedTetrahedronMass[tetrahedronIndex][index];
				}
			}
		} else if (location==HighOrderTetrahedronSetTopologyContainer::EDGE) {
			const TetrahedraAroundEdge &tae=adaptiveContainer->getTetrahedraAroundEdge(elementIndex);
			for (l=0;l<tae.size();++l) {
				tetrahedronIndex=tae[l];
				// finds the mass associated with that point in the tetrahedron
				const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=
					adaptiveContainer->getTetrahedronDOFArray(tetrahedronIndex);									
				for (j=0;indexArray[j]!=i;++j);
				assert(j<indexArray.size());
				// now sum the mass for each row element of the  local mass matrix
				for (k=0;k<indexArray.size();++k) {
					if (k<j) {
						index=k*indexArray.size()-k*(k-1)/2+j-k;
					} else {
						index=j*indexArray.size()-j*(j-1)/2+k-j;
					}
					mass+=reducedTetrahedronMass[tetrahedronIndex][index];
				}
			}
		} else if (location==HighOrderTetrahedronSetTopologyContainer::EDGE) {
			const TetrahedraAroundEdge &tae=adaptiveContainer->getTetrahedraAroundEdge(elementIndex);
			for (l=0;l<tae.size();++l) {
				tetrahedronIndex=tae[l];
				// finds the mass associated with that point in the tetrahedron
				const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=
					adaptiveContainer->getTetrahedronDOFArray(tetrahedronIndex);
				for (j=0;indexArray[j]!=i;++j);
				assert(j<indexArray.size());
				// now sum the mass for each row element of the  local mass matrix
				for (k=0;k<indexArray.size();++k) {
					if (k<j) {
						index=k*indexArray.size()-k*(k-1)/2+j-k;
					} else {
						index=j*indexArray.size()-j*(j-1)/2+k-j;
					}
					mass+=reducedTetrahedronMass[tetrahedronIndex][index];
				}
			}
		} else if (location==HighOrderTetrahedronSetTopologyContainer::POINT) {
			const TetrahedraAroundVertex &tav=adaptiveContainer->getTetrahedraAroundVertex(elementIndex);
			for (l=0;l<tav.size();++l) {
				tetrahedronIndex=tav[l];
				// finds the mass associated with that point in the tetrahedron
				const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=
					adaptiveContainer->getTetrahedronDOFArray(tetrahedronIndex);
				for (j=0;indexArray[j]!=i;++j);
				assert(j<indexArray.size());
				// now sum the mass for each row element of the  local mass matrix
				for (k=0;k<indexArray.size();++k) {
					if (k<j) {
						index=k*indexArray.size()-k*(k-1)/2+j-k;
					} else {
						index=j*indexArray.size()-j*(j-1)/2+k-j;
					}
					mass+=reducedTetrahedronMass[tetrahedronIndex][index];
				}
			}
		}
		// allows an error of 2% between the real lumped matrix and the only currently stored
		error=(mass-vertexMass[i])/vertexMass[i];
		assert(fabs(error)<2e-2);
		totalMass+=vertexMass[i];
	}
	std::cerr<<"totalMass= "<<totalMass<<std::endl;
#endif
}
template <class DataTypes, class MassType>
void AdaptiveBezierHighOrderMeshMatrixMass<DataTypes, MassType>::handleTopologyChange(core::topology::Topology *topo)
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
						
						size_t nbDofs;
                        HighOrderDegreeType degree=this->highOrderTetraGeo->getTopologyContainer()->getDegree();
						size_t nbControlPoints=(degree+1)*(degree+2)*(degree+3)/6;
                        const MassVectorVector& tetrahedronMass = this->tetrahedronMassInfo.getValue();
						MassVectorVector& reducedTetrahedronMass = *reducedTetrahedronMassInfo.beginEdit();
						const sofa::helper::vector<unsigned int> &tetraModified=btdc->getTetrahedronArray();
						size_t j,k,l,m,n,index;
						
						AdaptiveHighOrderTetrahedronSetTopologyContainer::WeightedDOFArray wda[2];
						AdaptiveHighOrderTetrahedronSetTopologyContainer::WeightedDOF wd[2];
						// the set of vertices for which it is necessary to recomputed the lumped mass
						std::set<size_t> vertexModifiedSet;
	
						for (size_t i=0;i<tetraModified.size();++i) {
						
							// updates the reduced Mass Matrix by splitting the initial mass matrix into the active DOFs
							const MassVector &tetinfo=tetrahedronMass[tetraModified[i]];
							MassVector &reducedMass=reducedTetrahedronMass[tetraModified[i]];
							const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=
								adaptiveContainer->getTetrahedronDOFArray(tetraModified[i]);
							const AdaptiveHighOrderTetrahedronSetTopologyContainer::SeqWeightedDOFArray &swda=
								adaptiveContainer->getWeightedDOFArray(tetraModified[i]);
							nbDofs=indexArray.size();
							for (j=0;j<nbDofs;++j) {
								vertexModifiedSet.insert(indexArray[j]);
							}
							// resize and fill with zero entries
							reducedMass.resize(nbDofs*(nbDofs+1)/2);
							std::fill(reducedMass.begin(),
								reducedMass.end(),MassType(0));
							for ( l=0,j=0;j<nbControlPoints;++j) {
								wda[0]=swda[j];
#ifndef NDEBUG
								Real weight;
								// verify that the  weights sum to 1
								weight=(Real)0.0f;
								for (m=0;m<wda[0].size();++m) {
									weight+=wda[0][m].second;
								}
								assert(fabs(weight-1.0f)<1e-4);
#endif
								for ( k=j;k<nbControlPoints;++k,++l) {
									wda[1]=swda[k];
									// mass matrix to be split among the mapped DOF
									const MassType &mat=tetinfo[l];

	
									for (m=0;m<wda[0].size();++m) {
										for (n=0;n<wda[1].size();++n) {
											// sort the indides of the 2 weighted array
											if (wda[0][m].first<wda[1][n].first) {
												wd[0]=wda[0][m];
												wd[1]=wda[1][n];
											} else {
												wd[1]=wda[0][m];
												wd[0]=wda[1][n];
											}

											// formula to get the index of the offdiagonal term in the array of stiffnesses
											index=wd[0].first*nbDofs-wd[0].first*(wd[0].first-1)/2+(wd[1].first-wd[0].first);
											if((j!=k) && (wd[0].first==wd[1].first))
												reducedMass[index]+=2*mat*wd[0].second*wd[1].second; 
											else if ((j==k) &&  (wd[0].first!=wd[1].first))
												reducedMass[index]+=0.5*mat*wd[0].second*wd[1].second; 
											else
												reducedMass[index]+=mat*wd[0].second*wd[1].second; 
										
		
										


										}
									}

								}
							}
	
						}
						
						/// now update the lumped mass of some vertices
						/// sum the terms of the mass matrix for eaxh row
						std::set<size_t>::iterator itv;
						HighOrderTetrahedronSetTopologyContainer::HighOrderTetrahedronPointLocation location;
						size_t elementIndex,elementOffset,tetrahedronIndex;
                        helper::vector<MassType> &vertexMass = *(this->vertexMassInfo.beginEdit());
						MassType mass;
						for (itv=vertexModifiedSet.begin();itv!=vertexModifiedSet.end();++itv) {
							mass=(MassType)0.0f;
							adaptiveContainer->getLocationFromGlobalIndex(*itv,location,elementIndex,elementOffset);
							if (location==HighOrderTetrahedronSetTopologyContainer::TETRAHEDRON) {
								// finds the mass associated with that point in the tetrahedron
								const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=
										adaptiveContainer->getTetrahedronDOFArray(elementIndex);

								for (j=0;indexArray[j]!=(*itv);++j);
								assert(j<indexArray.size());
								// now sum the mass for each row element of the  local mass matrix
								for (k=0;k<indexArray.size();++k) {
									if (k<j) {
										index=k*indexArray.size()-k*(k-1)/2+j-k;

									} else {
										index=j*indexArray.size()-j*(j-1)/2+k-j;
									}
									mass+=reducedTetrahedronMass[elementIndex][index];
								}
							} else if (location==HighOrderTetrahedronSetTopologyContainer::TRIANGLE) {
								const TetrahedraAroundTriangle &tat=adaptiveContainer->getTetrahedraAroundTriangle(elementIndex);
								for (l=0;l<tat.size();++l) {
									tetrahedronIndex=tat[l];
									// finds the mass associated with that point in the tetrahedron
									const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=
										adaptiveContainer->getTetrahedronDOFArray(tetrahedronIndex);
									for (j=0;indexArray[j]!=(*itv);++j);
									assert(j<indexArray.size());
									// now sum the mass for each row element of the  local mass matrix
									for (k=0;k<indexArray.size();++k) {
										if (k<j) {
											index=k*indexArray.size()-k*(k-1)/2+j-k;
										} else {
											index=j*indexArray.size()-j*(j-1)/2+k-j;
										}
										mass+=reducedTetrahedronMass[tetrahedronIndex][index];
									}
								}
							} else if (location==HighOrderTetrahedronSetTopologyContainer::EDGE) {
								const TetrahedraAroundEdge &tae=adaptiveContainer->getTetrahedraAroundEdge(elementIndex);
								for (l=0;l<tae.size();++l) {
									tetrahedronIndex=tae[l];
									// finds the mass associated with that point in the tetrahedron
									const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=
										adaptiveContainer->getTetrahedronDOFArray(tetrahedronIndex);
									for (j=0;indexArray[j]!=(*itv);++j);
									assert(j<indexArray.size());
									// now sum the mass for each row element of the  local mass matrix
									for (k=0;k<indexArray.size();++k) {
										if (k<j) {
											index=k*indexArray.size()-k*(k-1)/2+j-k;
										} else {
											index=j*indexArray.size()-j*(j-1)/2+k-j;
										}
										mass+=reducedTetrahedronMass[tetrahedronIndex][index];
									}
								}
							} else if (location==HighOrderTetrahedronSetTopologyContainer::EDGE) {
								const TetrahedraAroundEdge &tae=adaptiveContainer->getTetrahedraAroundEdge(elementIndex);
								for (l=0;l<tae.size();++l) {
									tetrahedronIndex=tae[l];
									// finds the mass associated with that point in the tetrahedron
									const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=
										adaptiveContainer->getTetrahedronDOFArray(tetrahedronIndex);
									for (j=0;indexArray[j]!=(*itv);++j);
									assert(j<indexArray.size());
									// now sum the mass for each row element of the  local mass matrix
									for (k=0;k<indexArray.size();++k) {
										if (k<j) {
											index=k*indexArray.size()-k*(k-1)/2+j-k;
										} else {
											index=j*indexArray.size()-j*(j-1)/2+k-j;
										}
										mass+=reducedTetrahedronMass[tetrahedronIndex][index];
									}
								}
							} else if (location==HighOrderTetrahedronSetTopologyContainer::POINT) {
								const TetrahedraAroundVertex &tav=adaptiveContainer->getTetrahedraAroundVertex(elementIndex);
								for (l=0;l<tav.size();++l) {
									tetrahedronIndex=tav[l];
									// finds the mass associated with that point in the tetrahedron
									const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=
										adaptiveContainer->getTetrahedronDOFArray(tetrahedronIndex);
									for (j=0;indexArray[j]!=(*itv);++j);
									assert(j<indexArray.size());
									// now sum the mass for each row element of the  local mass matrix
									for (k=0;k<indexArray.size();++k) {
										if (k<j) {
											index=k*indexArray.size()-k*(k-1)/2+j-k;
										} else {
											index=j*indexArray.size()-j*(j-1)/2+k-j;
										}
										mass+=reducedTetrahedronMass[tetrahedronIndex][index];
									}
								}
							}
							vertexMass[*itv]=mass;
						}
                        this->vertexMassInfo.endEdit();
						reducedTetrahedronMassInfo.endEdit();
#ifndef NDEBUG
						checkConsistencyOfLumpedMass();
#endif
					}
			}
			++itBegin;
		}

	}
}


template <class DataTypes, class MassType>
void AdaptiveBezierHighOrderMeshMatrixMass<DataTypes, MassType>::clear()
{

	MassVectorVector& tetrahedronMass = *reducedTetrahedronMassInfo.beginEdit();

	tetrahedronMass.clear();

	reducedTetrahedronMassInfo.endEdit();
}







template <class DataTypes, class MassType>
 const typename  AdaptiveBezierHighOrderMeshMatrixMass<DataTypes, MassType>::MassVector &
	 AdaptiveBezierHighOrderMeshMatrixMass<DataTypes, MassType>::getBezierTetrahedronMassVector(const size_t i) const {
		 return reducedTetrahedronMassInfo.getValue()[i];
 }








} // namespace mass

} // namespace component

} // namespace sofa

#endif
