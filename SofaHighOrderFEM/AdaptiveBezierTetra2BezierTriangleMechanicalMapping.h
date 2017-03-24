#ifndef SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERTETRA2BEZIERTRIANGLEMECHANICALMAPPING_H
#define SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERTETRA2BEZIERTRIANGLEMECHANICALMAPPING_H

#include <sofa/core/Mapping.h>
#include "initHighOrderFEM.h"
#include <SofaBaseTopology/TopologyData.h>

namespace sofa { namespace core { namespace topology { class BaseMeshTopology; } } }
namespace sofa { namespace component { namespace topology { class AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping; } } }
namespace sofa { namespace component { namespace topology { template<typename  D>  class BezierTriangleSetGeometryAlgorithms; } } }
namespace sofa { namespace component { namespace topology {  class HighOrderTriangleSetTopologyContainer; } } }
namespace sofa
{

namespace component
{

namespace mapping
{

template <class TIn, class TOut>
class SOFA_HIGHORDER_FEM_API AdaptiveBezierTetra2BezierTriangleMechanicalMapping : public core::Mapping<TIn, TOut>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(AdaptiveBezierTetra2BezierTriangleMechanicalMapping,TIn,TOut), SOFA_TEMPLATE2(core::Mapping,TIn,TOut));

    typedef core::Mapping<TIn, TOut> Inherit;
    typedef TIn In;
    typedef TOut Out;


    typedef typename Out::VecCoord OutVecCoord;
    typedef typename Out::VecDeriv OutVecDeriv;
    typedef typename Out::Coord OutCoord;
    typedef typename Out::Deriv OutDeriv;
    typedef typename Out::MatrixDeriv OutMatrixDeriv;


    typedef typename In::VecCoord InVecCoord;
    typedef typename In::VecDeriv InVecDeriv;
    typedef typename In::Coord InCoord;
    typedef typename In::Deriv InDeriv;
    typedef typename In::MatrixDeriv InMatrixDeriv;
    typedef typename InCoord::value_type Real;
	typedef sofa::defaulttype::Vec<3,Real> TriangleBarycentricCoordinates;

	typedef sofa::helper::vector<size_t> LocalVertexArray;
	typedef sofa::helper::vector<LocalVertexArray> SeqLocalVertexArray;


protected:
    AdaptiveBezierTetra2BezierTriangleMechanicalMapping(core::State<In>* from = NULL, core::State<Out>* to = NULL);

    virtual ~AdaptiveBezierTetra2BezierTriangleMechanicalMapping();

	sofa::component::topology::PointData<sofa::helper::vector<OutCoord > > offsetPositionData;

	 class PointADT2BTHandler : public   sofa::component::topology::TopologyDataHandler< sofa::core::topology::Topology::Point, helper::vector<OutCoord> >
	 {
	 public:
	
		 PointADT2BTHandler(AdaptiveBezierTetra2BezierTriangleMechanicalMapping<TIn,TOut>* _mm,  sofa::component::topology::PointData<sofa::helper::vector<OutCoord > >* _data)
			 : sofa::component::topology::TopologyDataHandler<sofa::core::topology::Topology::Point, sofa::helper::vector<OutCoord> >(_data), mm(_mm) {}
		     /// Move a list of points
		 virtual void move( const sofa::helper::vector<unsigned int> &indexList,
			 const sofa::helper::vector< sofa::helper::vector< unsigned int > >& ancestors,
			 const sofa::helper::vector< sofa::helper::vector< double > >& coefs);

	 protected:
		 AdaptiveBezierTetra2BezierTriangleMechanicalMapping<TIn,TOut>* mm;
	 };
	 PointADT2BTHandler *pointHandler;
public:

    void init();

    void apply(const core::MechanicalParams *mparams, Data<OutVecCoord>& out, const Data<InVecCoord>& in);

    void applyJ(const core::MechanicalParams *mparams, Data<OutVecDeriv>& out, const Data<InVecDeriv>& in);

    void applyJT(const core::MechanicalParams *mparams, Data<InVecDeriv>& out, const Data<OutVecDeriv>& in);

    void applyJT(const core::ConstraintParams *cparams, Data<InMatrixDeriv>& out, const Data<OutMatrixDeriv>& in);

	/// Function handling all the events (if listening=true)
    virtual void handleEvent(sofa::core::objectmodel::Event* event);

	/** \brief Translates the TopologyChange objects from the source to the target.
	*
	* Translates each of the TopologyChange objects waiting in the source list so that they have a meaning and
	* reflect the effects of the first topology changes on the second topology.
	*
	*/
	virtual void updateTopologicalMappingTopDown();

protected:
	// the associated topological map 
	topology::AdaptiveHighOrderTetra2HighOrderTriangleTopologicalMapping* topoMap;
	// the output bezier triangle geometry algorithm object
	topology::HighOrderTriangleSetTopologyContainer *btstc;
	// currently used bezier degree for the input Bezier triangulation or tetrahedral mesh
	size_t degree;
	// array of rest Position used to interpolate displacement
	sofa::helper::vector<InCoord> restPositionArray;
	sofa::helper::vector<InCoord> offsetPositionArray;
	sofa::helper::vector<size_t>  vertexTetra2TrianArray;


};





} // namespace mapping

} // namespace component

} // namespace sofa

#endif
