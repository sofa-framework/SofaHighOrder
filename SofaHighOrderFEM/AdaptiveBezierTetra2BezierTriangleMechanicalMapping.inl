#ifndef SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERTETRA2BEZIERTRIANGLEMECHANICALMAPPING_INL
#define SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERTETRA2BEZIERTRIANGLEMECHANICALMAPPING_INL

#include <sofa/core/MultiVecId.h>
#include "AdaptiveBezierTetra2BezierTriangleMechanicalMapping.h"
#include <sofa/simulation/AnimateEndEvent.h>
#include "AdaptiveBezierTetra2BezierTriangleTopologicalMapping.h"
#include "AdaptiveBezierTetrahedronSetTopologyContainer.h"
#include "HighOrderTriangleSetTopologyContainer.h"
#include <SofaBaseTopology/CommonAlgorithms.h>
#include <SofaBaseTopology/TopologyData.inl>

namespace sofa
{

namespace component
{

namespace mapping
{

typedef topology::AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping::WeightedDOFArray WeightedDOFArray;
typedef sofa::component::topology::HighOrderTriangleSetTopologyContainer HighOrderTriangleSetTopologyContainer;
typedef sofa::component::topology::HighOrderTetrahedronSetTopologyContainer HighOrderTetrahedronSetTopologyContainer;
typedef topology::AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping::Tetrahedron Tetrahedron;
typedef topology::AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping::Vec4 Vec4;

template <class TIn, class TOut>
AdaptiveBezierTetra2BezierTriangleMechanicalMapping<TIn, TOut>::AdaptiveBezierTetra2BezierTriangleMechanicalMapping(core::State<In>* from, core::State<Out>* to)
    : Inherit(from, to)
	, offsetPositionData(initData(&offsetPositionData, "offsetPosition", "offset Position Array"))
	    , topoMap(NULL)
    , degree(0)
	,pointHandler(NULL)
{
	pointHandler= new PointADT2BTHandler(this,&offsetPositionData);
}

template <class TIn, class TOut>
AdaptiveBezierTetra2BezierTriangleMechanicalMapping<TIn, TOut>::~AdaptiveBezierTetra2BezierTriangleMechanicalMapping()
{
	delete (pointHandler);
}


template <class TIn, class TOut>
void AdaptiveBezierTetra2BezierTriangleMechanicalMapping<TIn, TOut>::init()
{

	this->getContext()->get(topoMap);
	if (!topoMap) {
		serr << "Could not find any AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping object"<<sendl;
		return;
	}
	this->toModel->getContext()->get(btstc);
	if (!btstc){
		topoMap->fromModel->getContext()->get(btstc);
		if (!btstc){
			serr << "Could not find any HighOrderTriangleSetTopologyContainer object"<<sendl;
			return;
		}
	}
	
	degree=btstc->getDegree();

	// resize bezierTesselationWeightArray
	this->toModel->resize(topoMap->nbDofs);
	this->Inherit::init();
    this->f_listening.setValue(true);
	// store rest position
	size_t i,j,k;


	// store the rest position from the tetrahedron mechanical state.
    const InVecCoord &restPosition = this->fromModel->read(core::ConstVecCoordId::restPosition())->getValue();
	for (i=0;i<btstc->getNumberOfTriangularPoints();++i) {
		restPositionArray.push_back(restPosition[topoMap->tetraWeightedDOFArray[i][0].first]);
	}
	size_t rank=btstc->getNumberOfTriangularPoints();
	if (degree>1) {
		for (i=0;i<btstc->getNumberOfEdges();++i) {
			for (j=0;j<(degree-1);++j,++rank) {
				WeightedDOFArray &wda=topoMap->tetraWeightedDOFArray[rank];
				assert (wda.size()==1);
				restPositionArray.push_back(restPosition[wda[0].first]);
			}
		}
	}
	if (degree>2) {
		for (i=0;i<btstc->getNumberOfTriangles();++i) {
			for (j=1;j<(size_t)(degree-1);++j) {
				for (k=1;k<(degree-j);++k,++rank) {
					WeightedDOFArray &wda=topoMap->tetraWeightedDOFArray[rank];
					assert(wda.size()==1);
					// the edge is full degree : just copy the position
					restPositionArray.push_back(restPosition[wda[0].first]);
					
				}
			}
		}
	}
	assert(restPositionArray.size()==(topoMap->nbDofs));
	// initialize offset position to 
	sofa::helper::vector<OutCoord >	& offsetPos = *(offsetPositionData.beginEdit());
	offsetPos.resize(topoMap->nbDofs);
	offsetPositionData.endEdit();
	// store in an array the correspondance between tetrahedron vertices indices and triangle vertex indices (essentially faster).
	vertexTetra2TrianArray.resize(topoMap->from_btstc->getNbPoints());
	std::map<size_t,size_t>::iterator itv;
	for (itv=topoMap->vertexTetra2TrianMap.begin();
		itv!=topoMap->vertexTetra2TrianMap.end();++itv) {
			vertexTetra2TrianArray[(*itv).first]=(*itv).second;
	}
	// register PointData

	sofa::component::topology::HighOrderTriangleSetTopologyContainer *to_btstc;
    this->getContext()->get(to_btstc);
	assert(to_btstc!=NULL);
	offsetPositionData.createTopologicalEngine(to_btstc,pointHandler);
	offsetPositionData.registerTopologicalData();

	
}

template <class TIn, class TOut>
void AdaptiveBezierTetra2BezierTriangleMechanicalMapping<TIn, TOut>::PointADT2BTHandler::move( const sofa::helper::vector<unsigned int> &indexList,
			 const sofa::helper::vector< sofa::helper::vector< unsigned int > >& ancestors,
			 const sofa::helper::vector< sofa::helper::vector< double > >& coefs) 
{
	size_t i,elementIndex,elementOffset;
	unsigned int degree=mm->btstc->getDegree();
	HighOrderTriangleSetTopologyContainer::HighOrderTrianglePointLocation location;

	if (mm->topoMap->d_useSurfaceExtrapolation.getValue()) 
	{
        const InVecCoord &tetraPosition = mm->fromModel->read(core::ConstVecCoordId::position())->getValue();
        const OutVecCoord &trianPosition = mm->toModel->read(core::ConstVecCoordId::position())->getValue();


		sofa::helper::vector<OutCoord >	& offsetPos = *(mm->offsetPositionData.beginEdit());
		for (i=0;i<indexList.size();++i) {
			mm->btstc->getLocationFromGlobalIndex(indexList[i], location,elementIndex, elementOffset) ;
			if (location==HighOrderTriangleSetTopologyContainer::EDGE) {
				sofa::core::topology::Edge e=mm->btstc->getEdge(elementIndex);
				assert(ancestors[i].size()==1);
				Real w=(Real)(elementOffset+1)/degree;
				OutCoord pos;
				pos=tetraPosition[ancestors[i][0]];
				pos-=mm->restPositionArray[indexList[i]];
				assert(pos.norm()<1e-10);
				offsetPos[indexList[i]]=tetraPosition[ancestors[i][0]]-mm->restPositionArray[indexList[i]]-			
					(trianPosition[e[0]]-mm->restPositionArray[e[0]])*(1-w)-
					(trianPosition[e[1]]-mm->restPositionArray[e[1]])*w;
				assert(offsetPos[indexList[i]].norm()<1e-6);
			} else if  (location==HighOrderTriangleSetTopologyContainer::TRIANGLE) {
				sofa::core::topology::Triangle tr=mm->btstc->getTriangle(elementIndex);
				assert(ancestors[i].size()==1);
				topology::TriangleIndexVector tbc;
				mm->btstc->getTriangleIndexVectorFromTriangleOffset(elementOffset,tbc);
				OutCoord pos;
				pos=tetraPosition[ancestors[i][0]];
				pos-=mm->restPositionArray[indexList[i]];
				assert(pos.norm()<1e-10);
				offsetPos[indexList[i]]=tetraPosition[ancestors[i][0]]-mm->restPositionArray[indexList[i]]-			
					(trianPosition[tr[0]]-mm->restPositionArray[tr[0]])*tbc[0]/(Real)degree-
					(trianPosition[tr[1]]-mm->restPositionArray[tr[1]])*tbc[1]/(Real)degree-
					(trianPosition[tr[2]]-mm->restPositionArray[tr[2]])*tbc[2]/(Real)degree;
				assert(offsetPos[indexList[i]].norm()<1e-6);
			}

		}
		mm->offsetPositionData.endEdit();
	}
}

template <class TIn, class TOut>
void AdaptiveBezierTetra2BezierTriangleMechanicalMapping<TIn, TOut>::apply(const core::MechanicalParams * /*mparams*/, Data<OutVecCoord>& dOut, const Data<InVecCoord>& dIn)
{
	helper::WriteAccessor< Data<OutVecCoord> > out = dOut;
	helper::ReadAccessor< Data<InVecCoord> > in = dIn;

	const sofa::helper::vector<OutCoord >	& offsetPos = offsetPositionData.getValue();
	bool use_surfaceExtrapolation=topoMap->d_useSurfaceExtrapolation.getValue();


	size_t i;
	// copy the points of underlying triangulation
	for (i=0;i<btstc->getNumberOfTriangularPoints();++i) {
		out[i]=in[topoMap->tetraWeightedDOFArray[i][0].first];
	}
	
	size_t rank=btstc->getNumberOfTriangularPoints();
	// copy the points on  the edges pf the Bezier patches
	if (degree>1) {
		size_t j,k;
		for (i=0;i<btstc->getNumberOfEdges();++i) {
			for (j=0;j<(degree-1);++j,++rank) {
				WeightedDOFArray &wda=topoMap->tetraWeightedDOFArray[rank];
				if (wda.size()==1) {
					// the edge is full degree : just copy the position
					out[rank]=in[wda[0].first];
				} else {
					// the edge is linear : for now only simply interpolate position
					assert(wda.size()==2);
			//		out[rank]=in[wda[0].first]*wda[0].second+in[wda[1].first]*wda[1].second;
					//use weights coefficients
					if (use_surfaceExtrapolation) {
						out[rank]=restPositionArray[rank]+offsetPos[rank]+
							(in[wda[0].first]-restPositionArray[vertexTetra2TrianArray[wda[0].first]])*wda[0].second+
							(in[wda[1].first]-restPositionArray[vertexTetra2TrianArray[wda[1].first]])*wda[1].second;
					} else {
						// vertex position is computed from the barycentrical coordinates from 2 tetrahedra
						const std::vector<Vec4> &coefVal=topoMap->controlPointsBarycentricCoord[rank];
						assert(coefVal.size()==2);
						std::pair<Tetrahedron,Tetrahedron>  tetPair=topoMap->surfaceEdgeIndexToTetrahedronPairMap[i];
						out[rank]=OutCoord();
						for (k=0;k<4;++k) {
							out[rank]+=in[tetPair.first[k]]*coefVal[0][k];
						}
						for (k=0;k<4;++k) {
							out[rank]+=in[tetPair.second[k]]*coefVal[1][k];
						}
					}

				}
			}
		}
	}
	if (degree>2) {
		size_t j,k,l;
		for (i=0;i<btstc->getNumberOfTriangles();++i) {
			for (j=1;j<(size_t)(degree-1);++j) {
				for (k=1;k<(degree-j);++k,++rank) {
					WeightedDOFArray &wda=topoMap->tetraWeightedDOFArray[rank];
					if (wda.size()==1) {
						// the edge is full degree : just copy the position
						out[rank]=in[wda[0].first];
					} else {
						// the triangle  is linear : for now only simply interpolate position
						assert(wda.size()==3);
						if (use_surfaceExtrapolation) {
							out[rank]=restPositionArray[rank]+offsetPos[rank]+
								(in[wda[0].first]-restPositionArray[vertexTetra2TrianArray[wda[0].first]])*wda[0].second+
								(in[wda[1].first]-restPositionArray[vertexTetra2TrianArray[wda[1].first]])*wda[1].second+
								(in[wda[2].first]-restPositionArray[vertexTetra2TrianArray[wda[2].first]])*wda[2].second;
						} else {
							// vertex position is computed from the barycentrical coordinates from 1 tetrahedron
							const std::vector<Vec4> &coefVal=topoMap->controlPointsBarycentricCoord[rank];
							assert(coefVal.size()==1);
							Tetrahedron tet=topoMap->surfaceTriangleIndexToTetrahedronMap[i];
							out[rank]=OutCoord();
							for (l=0;l<4;++l) {
								out[rank]+=in[tet[l]]*coefVal[0][l];
							}
							
						}

//						out[rank]=in[wda[0].first]*wda[0].second+in[wda[1].first]*wda[1].second+
//							in[wda[2].first]*wda[2].second;
					}
				}
			}
		}
	}
	

}

template <class TIn, class TOut>
void AdaptiveBezierTetra2BezierTriangleMechanicalMapping<TIn, TOut>::applyJ(const core::MechanicalParams * /*mparams*/, Data<OutVecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{
    if (!topoMap) return;
	
    helper::WriteAccessor< Data<OutVecDeriv> > out = dOut;
	helper::ReadAccessor< Data<InVecDeriv> > in = dIn;
	bool use_surfaceExtrapolation=topoMap->d_useSurfaceExtrapolation.getValue();

	
	size_t i;
	// copy the points of underlying triangulation
	for (i=0;i<btstc->getNumberOfTriangularPoints();++i) {
		out[i]=in[topoMap->tetraWeightedDOFArray[i][0].first];
	}
	
	size_t rank=btstc->getNumberOfTriangularPoints();
	// copy the points on  the edges pf the Bezier patches
	if (degree>1) {
		size_t j,k;
		for (i=0;i<btstc->getNumberOfEdges();++i) {
			for (j=0;j<(degree-1);++j,++rank) {
				WeightedDOFArray &wda=topoMap->tetraWeightedDOFArray[rank];
				if (wda.size()==1) {
					// the edge is full degree : just copy the position
					out[rank]=in[wda[0].first];
				} else {
					// the edge is linear : for now only simply interpolate position
					assert(wda.size()==2);
					if (use_surfaceExtrapolation) {
						out[rank]=in[wda[0].first]*wda[0].second+in[wda[1].first]*wda[1].second;
					} else {
						// vertex position is computed from the barycentrical coordinates from 2 tetrahedra
						const std::vector<Vec4> &coefVal=topoMap->controlPointsBarycentricCoord[rank];
						assert(coefVal.size()==2);
						std::pair<Tetrahedron,Tetrahedron>  tetPair=topoMap->surfaceEdgeIndexToTetrahedronPairMap[i];
						out[rank]=OutDeriv();
						for (k=0;k<4;++k) {
							out[rank]+=in[tetPair.first[k]]*coefVal[0][k];
						}
						for (k=0;k<4;++k) {
							out[rank]+=in[tetPair.second[k]]*coefVal[1][k];
						}
					}
				}
			}
		}
	}
	if (degree>2) {
		size_t j,k,l;
		for (i=0;i<btstc->getNumberOfTriangles();++i) {
			for (j=1;j<(size_t)(degree-1);++j) {
				for (k=1;k<(degree-j);++k,++rank) {
					WeightedDOFArray &wda=topoMap->tetraWeightedDOFArray[rank];
					if (wda.size()==1) {
						// the edge is full degree : just copy the position
						out[rank]=in[wda[0].first];
					} else {
						// the triangle  is linear : for now only simply interpolate position
						assert(wda.size()==3);
						if (use_surfaceExtrapolation) {
							out[rank]=in[wda[0].first]*wda[0].second+in[wda[1].first]*wda[1].second+
								in[wda[2].first]*wda[2].second;
						} else {
							// vertex position is computed from the barycentrical coordinates from 1 tetrahedron
							const std::vector<Vec4> &coefVal=topoMap->controlPointsBarycentricCoord[rank];
							assert(coefVal.size()==1);
							Tetrahedron tet=topoMap->surfaceTriangleIndexToTetrahedronMap[i];
							out[rank]=OutCoord();
							for (l=0;l<4;++l) {
								out[rank]+=in[tet[l]]*coefVal[0][l];
							}
						}
					}
				}
			}
		}
	}
	


}

template <class TIn, class TOut>
void AdaptiveBezierTetra2BezierTriangleMechanicalMapping<TIn, TOut>::applyJT(const core::MechanicalParams * /*mparams*/, Data<InVecDeriv>& dOut, const Data<OutVecDeriv>& dIn)
{
    if (!topoMap) return;

    helper::WriteAccessor< Data<InVecDeriv> > out = dOut;
	helper::ReadAccessor< Data<OutVecDeriv> > in = dIn;
	bool use_surfaceExtrapolation=topoMap->d_useSurfaceExtrapolation.getValue();
	
	
	size_t i;
	// copy the points of underlying triangulation
	for (i=0;i<btstc->getNumberOfTriangularPoints();++i) {
		out[topoMap->tetraWeightedDOFArray[i][0].first]+=in[i];
	}
	size_t rank=btstc->getNumberOfTriangularPoints();
	// copy the points on  the edges pf the Bezier patches
	if (degree>1) {
		size_t j,k;
		for (i=0;i<btstc->getNumberOfEdges();++i) {
			for (j=0;j<(degree-1);++j,++rank) {
				WeightedDOFArray &wda=topoMap->tetraWeightedDOFArray[rank];
				if (wda.size()==1) {
					// the edge is full degree : just copy the position
					out[wda[0].first]+=in[rank];
				} else {
					// the edge is linear : for now only simply interpolate position
					assert(wda.size()==2);
					if (use_surfaceExtrapolation) {
						out[wda[0].first]+=in[rank]*wda[0].second;
						out[wda[1].first]+=in[rank]*wda[1].second;
					} else {
						// vertex position is computed from the barycentrical coordinates from 2 tetrahedra
						const std::vector<Vec4> &coefVal=topoMap->controlPointsBarycentricCoord[rank];
						assert(coefVal.size()==2);
						std::pair<Tetrahedron,Tetrahedron>  tetPair=topoMap->surfaceEdgeIndexToTetrahedronPairMap[i];
						//out[rank]=OutDeriv();
						for (k=0;k<4;++k) {
							out[tetPair.first[k]]+=in[rank]*coefVal[0][k];
						}
						for (k=0;k<4;++k) {
							out[tetPair.second[k]]+=in[rank]*coefVal[1][k];
						}
					}
				}
			}
		}
	}
	if (degree>2) {
		size_t j,k,l;
		for (i=0;i<btstc->getNumberOfTriangles();++i) {
			for (j=1;j<(size_t)(degree-1);++j) {
				for (k=1;k<(degree-j);++k,++rank) {
					WeightedDOFArray &wda=topoMap->tetraWeightedDOFArray[rank];
					if (wda.size()==1) {
						// the edge is full degree : just copy the position
						out[wda[0].first]+=in[rank];
					} else {
						// the triangle  is linear : for now only simply interpolate position
						assert(wda.size()==3);
						if (use_surfaceExtrapolation) {
							out[wda[0].first]+=in[rank]*wda[0].second;
							out[wda[1].first]+=in[rank]*wda[1].second;
							out[wda[2].first]+=in[rank]*wda[2].second;
						} else {
							// vertex position is computed from the barycentrical coordinates from 1 tetrahedron
							const std::vector<Vec4> &coefVal=topoMap->controlPointsBarycentricCoord[rank];
							assert(coefVal.size()==1);
							Tetrahedron tet=topoMap->surfaceTriangleIndexToTetrahedronMap[i];
							// out[rank]=OutCoord();
							for (l=0;l<4;++l) {
								out[tet[l]]+=in[rank]*coefVal[0][l];
							}
						}
					}
				}
			}
		}
	}
}


template <class TIn, class TOut>
void AdaptiveBezierTetra2BezierTriangleMechanicalMapping<TIn, TOut>::updateTopologicalMappingTopDown()
{

    // INITIALISATION of TRIANGULAR mesh from TETRAHEDRAL mesh :
//	cerr << "updateTopologicalMappingTopDown called" << endl;

    if (this->fromModel)
    {

     //   TriangleSetTopologyModifier *to_tstm;
    //    toModel->getContext()->get(to_tstm);

        if (this->toModel)
		{

			std::list<const sofa::core::topology::TopologyChange *>::const_iterator itBegin=topoMap->fromModel->beginChange();
			std::list<const  sofa::core::topology::TopologyChange *>::const_iterator itEnd=topoMap->fromModel->endChange();

			//sofa::helper::vector <unsigned int>& Loc2GlobVec = *(Loc2GlobDataVec.beginEdit());

			while( itBegin != itEnd )
			{
				 sofa::core::topology::TopologyChangeType changeType = (*itBegin)->getChangeType();

				switch( changeType )
				{

				case core::topology::POINTSREMOVED:
					{
						std::cerr<<" got points removed"<<std::endl;
						break;
					}
				}
			}
		}
	}
}
template <class TIn, class TOut>
void AdaptiveBezierTetra2BezierTriangleMechanicalMapping<TIn, TOut>::handleEvent(sofa::core::objectmodel::Event* event)
{
	if (dynamic_cast<sofa::simulation::AnimateEndEvent *>(event))
	{
	/*	const sofa::core::MechanicalParams* mparams=sofa::core::MechanicalParams::defaultInstance();
		sofa::core::State<TOut>*  toModel = this->toModel.get(mparams);
		sofa::core::State<TIn>* fromModel = this->fromModel.get(mparams);
		sofa::core::MultiVecCoordId outPos=sofa::core::VecCoordId::position();
		sofa::core::ConstMultiVecCoordId inPos= sofa::core::ConstVecCoordId::position();
		OutDataVecCoord* out = outPos[toModel].write();
		const InDataVecCoord* in = inPos[fromModel].read();
		this->apply(mparams,*out,*in );
		sofa::core::MultiVecDerivId outVel=sofa::core::VecDerivId::velocity();
		sofa::core::ConstMultiVecDerivId inVel=sofa::core:: ConstVecDerivId::velocity();
		 OutDataVecDeriv* outV = outVel[toModel].write();
        const InDataVecDeriv* inV = inVel[fromModel].read();
		this->applyJ(mparams,*outV,*inV ); */
//		applyJ(sofa::core::MechanicalParams::defaultInstance(), sofa::core::VecDerivId::velocity(), sofa::core::ConstVecDerivId::velocity());
	}
}

template <class TIn, class TOut>
void AdaptiveBezierTetra2BezierTriangleMechanicalMapping<TIn, TOut>::applyJT(const core::ConstraintParams * /*cparams*/, Data<InMatrixDeriv>& /*dOut*/, const Data<OutMatrixDeriv>& /*dIn*/)
{

    if (!topoMap)
        return;

//    const sofa::helper::vector< std::pair< Mesh2PointTopologicalMapping::Element, int> >& pointSource = topoMap->getPointSource();

  //  if (pointSource.empty())
   //     return;
	/*
    InMatrixDeriv& out = *dOut.beginEdit();
    const OutMatrixDeriv& in = dIn.getValue();

    const core::topology::BaseMeshTopology::SeqEdges& edges = inputTopo->getEdges();
    const core::topology::BaseMeshTopology::SeqTriangles& triangles = inputTopo->getTriangles();
    const core::topology::BaseMeshTopology::SeqQuads& quads = inputTopo->getQuads();
    const core::topology::BaseMeshTopology::SeqTetrahedra& tetrahedra = inputTopo->getTetrahedra();
    const core::topology::BaseMeshTopology::SeqHexahedra& hexahedra = inputTopo->getHexahedra();

    typename Out::MatrixDeriv::RowConstIterator rowItEnd = in.end();

    for (typename Out::MatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
    {
        typename Out::MatrixDeriv::ColConstIterator colIt = rowIt.begin();
        typename Out::MatrixDeriv::ColConstIterator colItEnd = rowIt.end();

        // Creates a constraints if the input constraint is not empty.
        if (colIt != colItEnd)
        {
            typename In::MatrixDeriv::RowIterator o = out.writeLine(rowIt.index());

            while (colIt != colItEnd)
            {
                const unsigned int indexIn = colIt.index();
                const OutDeriv data = colIt.val();
                std::pair< Mesh2PointTopologicalMapping::Element, int> source = pointSource[indexIn];

                switch (source.first)
                {
                case topology::Mesh2PointTopologicalMapping::POINT:
                {
                    o.addCol(source.second, data);

                    break;
                }
                case topology::Mesh2PointTopologicalMapping::EDGE:
                {
                    core::topology::BaseMeshTopology::Edge e = edges[source.second];
                    typename In::Deriv f = data;
                    double fx = topoMap->getEdgeBaryCoords()[indexIn][0];

                    o.addCol(e[0], f * (1 - fx));
                    o.addCol(e[1], f * fx);

                    break;
                }
                case topology::Mesh2PointTopologicalMapping::TRIANGLE:
                {
                    core::topology::BaseMeshTopology::Triangle t = triangles[source.second];
                    typename In::Deriv f = data;
                    double fx = topoMap->getTriangleBaryCoords()[indexIn][0];
                    double fy = topoMap->getTriangleBaryCoords()[indexIn][1];

                    o.addCol(t[0], f * (1 - fx - fy));
                    o.addCol(t[1], f * fx);
                    o.addCol(t[2], f * fy);

                    break;
                }
                case topology::Mesh2PointTopologicalMapping::QUAD:
                {
                    core::topology::BaseMeshTopology::Quad q = quads[source.second];
                    typename In::Deriv f = data;
                    double fx = topoMap->getQuadBaryCoords()[indexIn][0];
                    double fy = topoMap->getQuadBaryCoords()[indexIn][1];

                    o.addCol(q[0], f * ((1-fx) * (1-fy)));
                    o.addCol(q[1], f * (fx * (1-fy)));
					o.addCol(q[2], f * ((1-fx) * fy));
					o.addCol(q[3], f * (fx * fy));

                    break;
                }
                case topology::Mesh2PointTopologicalMapping::TETRA:
                {
                    core::topology::BaseMeshTopology::Tetra t = tetrahedra[source.second];
                    typename In::Deriv f = data;
                    double fx = topoMap->getTetraBaryCoords()[indexIn][0];
                    double fy = topoMap->getTetraBaryCoords()[indexIn][1];
                    double fz = topoMap->getTetraBaryCoords()[indexIn][2];

                    o.addCol(t[0], f * (1-fx-fy-fz));
                    o.addCol(t[1], f * fx);
					o.addCol(t[2], f * fy);
					o.addCol(t[3], f * fz);

                    break;
                }
                case topology::Mesh2PointTopologicalMapping::HEXA:
                {
                    core::topology::BaseMeshTopology::Hexa h = hexahedra[source.second];
                    typename In::Deriv f = data;
                    const double fx = topoMap->getHexaBaryCoords()[indexIn][0];
                    const double fy = topoMap->getHexaBaryCoords()[indexIn][1];
                    const double fz = topoMap->getHexaBaryCoords()[indexIn][2];
                    const double oneMinFx = 1 - fx;
                    const double oneMinFy = 1 - fy;
                    const double oneMinFz = 1 - fz;

                    o.addCol(h[0] , f * oneMinFx * oneMinFy * oneMinFz);
                    o.addCol(h[1] , f * fx * oneMinFy * oneMinFz);
					o.addCol(h[3] , f * oneMinFx * fy * oneMinFz);
					o.addCol(h[2] , f * fx * fy * oneMinFz);
					o.addCol(h[4] , f * oneMinFx * oneMinFy * fz);
					o.addCol(h[5] , f * fx * oneMinFy * fz);
					o.addCol(h[6] , f * fx * fy * fz);
					o.addCol(h[7] , f * oneMinFx * fy * fz);

                    break;
                }
                default:

                    break;
                }

                ++colIt;
            }
        }
    }
	
    dOut.endEdit(); */
}

} // namespace mapping

} // namespace component

} // namespace sofa

#endif
