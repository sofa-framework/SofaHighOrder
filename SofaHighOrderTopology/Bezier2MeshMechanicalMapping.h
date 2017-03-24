
#ifndef SOFA_HIGHORDERTOPOLOGY_BEZIER2MESHMECHANICALMAPPING_H
#define SOFA_HIGHORDERTOPOLOGY_BEZIER2MESHMECHANICALMAPPING_H

#include <sofa/core/Mapping.h>

#include <sofa/defaulttype/VecTypes.h>

namespace sofa { namespace core { namespace topology { class BaseMeshTopology; } } }
namespace sofa { namespace component { namespace topology { class HighOrder2MeshTopologicalMapping; } } }
namespace sofa { namespace component { namespace topology { template<typename  D>  class BezierTriangleSetGeometryAlgorithms; } } }
namespace sofa { namespace component { namespace topology {  class HighOrderTriangleSetTopologyContainer; } } }
namespace sofa
{

namespace component
{

namespace mapping
{

template <class TIn, class TOut>
class Bezier2MeshMechanicalMapping : public core::Mapping<TIn, TOut>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(Bezier2MeshMechanicalMapping,TIn,TOut), SOFA_TEMPLATE2(core::Mapping,TIn,TOut));

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
    typedef sofa::defaulttype::Vec<3,Real> Vec3;
protected:
    Bezier2MeshMechanicalMapping(core::State<In>* from = NULL, core::State<Out>* to = NULL);

    virtual ~Bezier2MeshMechanicalMapping();

public:

    void init();

    void apply(const core::MechanicalParams *mparams, Data<OutVecCoord>& out, const Data<InVecCoord>& in);

    void applyJ(const core::MechanicalParams *mparams, Data<OutVecDeriv>& out, const Data<InVecDeriv>& in);

    void applyJT(const core::MechanicalParams *mparams, Data<InVecDeriv>& out, const Data<OutVecDeriv>& in);

    void applyJT(const core::ConstraintParams *cparams, Data<InMatrixDeriv>& out, const Data<OutMatrixDeriv>& in);

protected:
	// the associated topological map 
	topology::HighOrder2MeshTopologicalMapping* topoMap;

	// the input bezier triangle geometry algorithm object
	topology::HighOrderTriangleSetTopologyContainer *btstc;
	// currently used bezier degree for the input Bezier triangulation or tetrahedral mesh
	size_t bezierDegree;
	// currently used tesselation degree for the output riangulation or tetrahedral mesh
	size_t tesselationDegree;

	/// precomputed coefficients to interpolate the positions of points.
	sofa::helper::vector< sofa::helper::vector<Real> > precomputedLinearBernsteinCoefficientArray; 
	sofa::helper::vector< sofa::helper::vector<Real> > precomputedTriangularBernsteinCoefficientArray; 
	sofa::helper::vector< sofa::helper::vector<Real> > precomputedDerivUTriangularBernsteinCoefficientArray; 
	sofa::helper::vector< sofa::helper::vector<Real> > precomputedDerivVTriangularBernsteinCoefficientArray;
	/// precompute weight array for tesselated mesh
	sofa::helper::vector< Real > bezierTesselationWeightArray; 
	// local indexing of points inside tessellated triangles
	sofa::helper::vector<sofa::defaulttype::Vec<3,unsigned char > > tesselatedTriangleIndices; 


};



#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_HIGHORDERTOPOLOGY_BEZIER2MESHMECHANICALMAPPING_CPP)  //// ATTENTION PB COMPIL WIN3Z
#ifndef SOFA_FLOAT
extern template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< defaulttype::Vec3dTypes, defaulttype::Vec3dTypes >;
extern template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< defaulttype::Vec1dTypes, defaulttype::Vec1dTypes >;
extern template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< defaulttype::Vec3dTypes, defaulttype::ExtVec3dTypes >;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< defaulttype::Vec3fTypes, defaulttype::Vec3fTypes >;
extern template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< defaulttype::Vec1fTypes, defaulttype::Vec1fTypes >;
extern template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< defaulttype::Vec3fTypes, defaulttype::ExtVec3fTypes >;
#endif

#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
extern template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< defaulttype::Vec3dTypes, defaulttype::Vec3fTypes >;
extern template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< defaulttype::Vec3fTypes, defaulttype::Vec3dTypes >;
extern template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< defaulttype::Vec1dTypes, defaulttype::Vec1dTypes >;
extern template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< defaulttype::Vec1fTypes, defaulttype::Vec1fTypes >;
extern template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< defaulttype::Vec3fTypes, defaulttype::ExtVec3dTypes >;
extern template class SOFA_HIGHORDER_TOPOLOGY_API Bezier2MeshMechanicalMapping< defaulttype::Vec3dTypes, defaulttype::ExtVec3fTypes >;
#endif
#endif
#endif

} // namespace mapping

} // namespace component

} // namespace sofa

#endif
