
#ifndef SOFA_HIGHORDERDFEM_HIGHORDERMESHMATRIXMASS_H
#define SOFA_HIGHORDERDFEM_HIGHORDERMESHMATRIXMASS_H
#include "initHighOrderFEM.h"

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/behavior/Mass.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/Event.h>
#include <SofaBaseTopology/TopologyData.h>
#include <sofa/helper/vector.h>
#include <sofa/defaulttype/RigidTypes.h>
//VERY IMPORTANT FOR GRAPHS
#include <sofa/helper/map.h>

#include <sofa/core/topology/BaseMeshTopology.h>

namespace sofa
{
namespace component
{
namespace topology
{
	/// forward declaration to avoid adding includes in .h
	template< class DataTypes> class EdgeSetGeometryAlgorithms;
	template< class DataTypes> class TriangleSetGeometryAlgorithms;
	template< class DataTypes> class TetrahedronSetGeometryAlgorithms;
	template< class DataTypes> class HighOrderTetrahedronSetGeometryAlgorithms;
	template< class DataTypes> class HighOrderTriangleSetGeometryAlgorithms;
	template< class DataTypes> class QuadSetGeometryAlgorithms;
	template< class DataTypes> class HexahedronSetGeometryAlgorithms;
}

namespace mass
{

template<class DataTypes, class TMassType>
class HighOrderMeshMatrixMassInternalData
{
};



// template<class Vec> void readVec1(Vec& vec, const char* str);
template <class DataTypes, class TMassType>
class SOFA_HIGHORDER_FEM_API HighOrderMeshMatrixMass : public core::behavior::Mass<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(HighOrderMeshMatrixMass,DataTypes,TMassType), SOFA_TEMPLATE(core::behavior::Mass,DataTypes));

    typedef core::behavior::Mass<DataTypes> Inherited;
    typedef typename DataTypes::VecCoord                    VecCoord;
    typedef typename DataTypes::VecDeriv                    VecDeriv;
    typedef typename DataTypes::Coord                       Coord;
    typedef typename DataTypes::Deriv                       Deriv;
    typedef typename DataTypes::Real                        Real;
    typedef core::objectmodel::Data<VecCoord>               DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv>               DataVecDeriv;
    typedef TMassType                                       MassType;
    typedef helper::vector<MassType>                        MassVector;
    typedef helper::vector<MassVector>                      MassVectorVector;
	typedef sofa::defaulttype::Vec<2,size_t>				TetrahedronIndexVectorIndexPair;
	typedef std::pair<size_t,Real>         RegularMassCoefficientEntry;
	typedef sofa::core::topology::Topology::Tetrahedron Tetrahedron;
	typedef sofa::core::topology::Topology::TetraID TetraID;
	typedef sofa::core::topology::Topology::Tetra Tetra;
	typedef sofa::core::topology::Topology::Point Point;
	typedef sofa::core::topology::Topology::Triangle Triangle;
	typedef sofa::core::topology::Topology::Edge Edge;
	typedef sofa::core::topology::Topology::Quad Quad;
	typedef sofa::core::topology::Topology::Hexahedron Hexahedron;
	typedef sofa::core::topology::BaseMeshTopology::EdgesInTriangle EdgesInTriangle;
	typedef sofa::core::topology::BaseMeshTopology::EdgesInTetrahedron EdgesInTetrahedron;
	typedef sofa::core::topology::BaseMeshTopology::EdgesInQuad EdgesInQuad;
	typedef sofa::core::topology::BaseMeshTopology::EdgesInHexahedron EdgesInHexahedron;
	typedef sofa::core::topology::BaseMeshTopology::TrianglesInTetrahedron TrianglesInTetrahedron;


    // In case of non 3D template
    typedef defaulttype::Vec<2,Real> Vec2;
    typedef defaulttype::Vec<3,Real> Vec3;
    typedef defaulttype::Vec<4,Real> Vec4;
    /// assumes the geometry object type is 3D
    typedef defaulttype::StdVectorTypes< Vec3, Vec3, Real > GeometricalTypes ;
    typedef defaulttype::StdVectorTypes< Vec2, Vec2, Real > GeometricalTypes2D ;

    /// Topological enum to classify encounter meshes
    typedef enum
    {
        TOPOLOGY_UNKNOWN=0,
        TOPOLOGY_EDGESET=1,
        TOPOLOGY_TRIANGLESET=2,
        TOPOLOGY_TETRAHEDRONSET=3,
        TOPOLOGY_QUADSET=4,
        TOPOLOGY_HEXAHEDRONSET=5,
        TOPOLOGY_BEZIERTETRAHEDRONSET=6,
        TOPOLOGY_BEZIERTRIANGLESET=7
    } TopologyType;
	/// the way the mass should be computed on non-linear elements
	typedef enum 
	{
		EXACT_INTEGRATION=1,
		NUMERICAL_INTEGRATION=2,
		AFFINE_ELEMENT_INTEGRATION=3,
        BEZIER_NUMERICAL_INTEGRATION = 4
	} IntegrationMethod;
	/// the way the mass should be lumped as a diagonal matrix
	typedef enum 
	{
		ROW_SUM=1,
		SCALED_DIAGONAL=2
	} LumpingMethod;

    /// Mass info are stocked on vertices and edges (if lumped matrix)
    topology::PointData<helper::vector<MassType> >  vertexMassInfo;
    topology::EdgeData<helper::vector<MassType> >   edgeMassInfo;

    /* ---------- Specific data for Bezier Elements ------*/
    /// use this data structure to store mass for Bezier tetrahedra. 
    //// The size of the vector is nbControlPoints*(nbControlPoints+1)/2 where nbControlPoints=(degree+1)*(degree+2)*(degree+3)/2
    topology::TetrahedronData<helper::vector<MassVector> > tetrahedronMassInfo;
    topology::TriangleData<helper::vector<MassVector> > triangleMassInfo;
    // array of Tetrahedral Bezier indices
    //sofa::helper::vector<TetrahedronIndexVector> tbiArray;
    /* ---------- end ------*/

    /// the mass density used to compute the mass from a mesh topology and geometry
    Data< Real >         m_massDensity;

    /// to display the center of gravity of the system
    Data< bool >         showCenterOfGravity;
    Data< Real >         showAxisSize;
    /// if mass lumping should be performed (only compute mass on vertices)
    Data< bool >         lumping;
    /// if specific mass information should be outputed
    Data< bool >         printMass;
    Data<std::map < std::string, sofa::helper::vector<double> > > f_graph;
    /// the order of integration for numerical integration
    Data<size_t>	     numericalIntegrationOrder;
    /// the type of numerical integration method chosen
    Data<std::string>	     numericalIntegrationMethod;
    /// the type of integration method chosen for non linear element.
    Data<std::string>	 d_integrationMethod; 

	// measure the time spent in assembling the mass matrix
	Data<Real> d_assemblyTime;
	// whether each affine element should be assembled with the affine assembly method irrespective to the chosen integration method  
	Data<bool> d_forceAffineAssemblyForAffineElements;
    /// the type of integration method chosen for non linear element.
    Data<std::string>	 d_lumpingMethod; 
    LumpingMethod    lumpingMethod;
    // the type of integration method for the mass matrix
    IntegrationMethod    integrationMethod;

protected:

    /// The type of topology to build the mass from the topology
    TopologyType topologyType;
    Real massLumpingCoeff;
    Real savedMass;

    HighOrderMeshMatrixMass();
    ~HighOrderMeshMatrixMass();

    /// Internal data required for Cuda computation (copy of vertex mass for deviceRead)
    HighOrderMeshMatrixMassInternalData<DataTypes, MassType> data;
    friend class HighOrderMeshMatrixMassInternalData<DataTypes, MassType>;

public:

    sofa::core::topology::BaseMeshTopology* _topology;

    sofa::component::topology::EdgeSetGeometryAlgorithms<GeometricalTypes>* edgeGeo;
    sofa::component::topology::TriangleSetGeometryAlgorithms<GeometricalTypes>* triangleGeo;
    sofa::component::topology::QuadSetGeometryAlgorithms<GeometricalTypes>* quadGeo;
    sofa::component::topology::TetrahedronSetGeometryAlgorithms<GeometricalTypes>* tetraGeo;
    sofa::component::topology::HexahedronSetGeometryAlgorithms<GeometricalTypes>* hexaGeo;
    sofa::component::topology::HighOrderTetrahedronSetGeometryAlgorithms<GeometricalTypes>* highOrderTetraGeo;
    sofa::component::topology::HighOrderTriangleSetGeometryAlgorithms<GeometricalTypes2D>* highOrderTrianGeo;

    virtual void clear();

    virtual void reinit();
    virtual void init();

    TopologyType getMassTopologyType() const
    {
        return topologyType;
    }

    void setMassTopologyType(TopologyType t)
    {
        topologyType = t;
    }


    Real getMassDensity() const
    {
        return m_massDensity.getValue();
    }

    void setMassDensity(Real m)
    {
        m_massDensity.setValue(m);
    }

    /// Copy the vertex mass scalar (in case of CudaTypes)
    void copyVertexMass();


    // -- Mass interface
    void addMDx(const core::MechanicalParams*, DataVecDeriv& f, const DataVecDeriv& dx, SReal factor);

    void accFromF(const core::MechanicalParams*, DataVecDeriv& a, const DataVecDeriv& f); // This function can't be used as it use M^-1

    void addForce(const core::MechanicalParams*, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v);

    SReal getKineticEnergy(const core::MechanicalParams*, const DataVecDeriv& v) const;  ///< vMv/2 using dof->getV()

    SReal getPotentialEnergy(const core::MechanicalParams*, const DataVecCoord& x) const;   ///< Mgx potential in a uniform gravity field, null at origin

    defaulttype::Vector6 getMomentum(const core::MechanicalParams* mparams, const DataVecCoord& x, const DataVecDeriv& v) const;  ///< (Mv,cross(x,Mv))

    void addGravityToV(const core::MechanicalParams* mparams, DataVecDeriv& d_v);

    bool isDiagonal() {return false;}



    /// Add Mass contribution to global Matrix assembling
    void addMToMatrix(const core::MechanicalParams *mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix);

    SReal getElementMass(unsigned int index) const;
    void getElementMass(unsigned int index, defaulttype::BaseMatrix *m) const;

    void draw(const core::visual::VisualParams* vparams);

    /// Answer wether mass matrix is lumped or not
    bool isLumped() { return lumping.getValue(); }
	// returns the mass vector for a given index of a Bezier tetrahedron
	virtual const  MassVector &getBezierTetrahedronMassVector(const size_t i) const;
	// returns the mass vector for a given index of a Bezier triangle
	virtual const  MassVector &getBezierTriangleMassVector(const size_t i) const;

protected:
	// the array where mass coefficients are stored for affine elements.
	std::vector<Real> affineMassCoefficientArray;

	// the array where regular mass coefficients are stored for regular elements.
	std::vector<Real> regularMassCoefficientArray;
	// the data stored for each integration point
	struct NumericalIntegrationMassData {

		// for each pair of control point store the weight matrix  w_\gamma * density * N_p (\param_gamma) *  N_q(\param_gamma)
		std::vector<Real> weightArray;
		// for each control point store the weight matrix  w_\gamma * density * N_p (\param_gamma) 
		std::vector<Real> weightLumpedArray;
		/// barycentric coordinate of the integration point \param_i
        // tetrahedral barycentric coordinates of integration point
		Vec4 integrationPoint;
        // triangular barycentric coordinates of integration point
		Vec3 integrationPointTriangle;
		Real weight;
	};
    // in the BEZIER_NUMERICAL_INTEGRATION store correspondance between pair of shape function of degree d and shape functions of degree 2d.
    std::vector<size_t> indexCorrespondance;
    // for each pair of control point store the weight equal to 1/(p!q!) 
    std::vector<Real> weightBezierArray;

    class VertexMassHandler : public topology::TopologyDataHandler<Point,MassVector>
    {
    public:
        VertexMassHandler(HighOrderMeshMatrixMass<DataTypes,TMassType>* _m, topology::PointData<helper::vector<TMassType> >* _data) : topology::TopologyDataHandler<Point,helper::vector<TMassType> >(_data), m(_m) {}

        /// Mass initialization Creation Functions:
        /// Vertex mass coefficient matrix creation function
        void applyCreateFunction(unsigned int pointIndex, TMassType & VertexMass,
                const sofa::helper::vector< unsigned int > &,
                const sofa::helper::vector< double >&);


        ///////////////////////// Functions on Triangles //////////////////////////////////////

        /// Mass coefficient Creation/Destruction functions for Triangular Mesh:
        /// Vertex coefficient of mass matrix creation function to handle creation of new triangles
        void applyTriangleCreation(const sofa::helper::vector< unsigned int >& /*indices*/,
                const sofa::helper::vector< Triangle >& /*elems*/,
                const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
                const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/);

        /// Vertex coefficient of mass matrix destruction function to handle creation of new triangles
        void applyTriangleDestruction(const sofa::helper::vector<unsigned int> & /*indices*/);

        using topology::TopologyDataHandler<Point,MassVector>::ApplyTopologyChange;
        /// Callback to add triangles elements.
        void ApplyTopologyChange(const core::topology::TrianglesAdded* /*event*/);
        /// Callback to remove triangles elements.
        void ApplyTopologyChange(const core::topology::TrianglesRemoved* /*event*/);


        ///////////////////////// Functions on Quads //////////////////////////////////////

        /// Mass coefficient Creation/Destruction functions for Quad Mesh:
        /// Vertex coefficient of mass matrix creation function to handle creation of new quads
        void applyQuadCreation(const sofa::helper::vector< unsigned int >& /*indices*/,
                const sofa::helper::vector< Quad >& /*elems*/,
                const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
                const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/);

        /// Vertex coefficient of mass matrix destruction function to handle creation of new quads
        void applyQuadDestruction(const sofa::helper::vector<unsigned int> & /*indices*/);

        /// Callback to add quads elements.
        void ApplyTopologyChange(const core::topology::QuadsAdded* /*event*/);
        /// Callback to remove quads elements.
        void ApplyTopologyChange(const core::topology::QuadsRemoved* /*event*/);


        ///////////////////////// Functions on Tetrahedron //////////////////////////////////////

        /// Mass coefficient Creation/Destruction functions for Tetrahedral Mesh:
        /// Vertex coefficient of mass matrix creation function to handle creation of new tetrahedra
        void applyTetrahedronCreation(const sofa::helper::vector< unsigned int >& /*indices*/,
                const sofa::helper::vector< Tetrahedron >& /*elems*/,
                const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
                const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/);

        /// Vertex coefficient of mass matrix destruction function to handle creation of new tetrahedra
        void applyTetrahedronDestruction(const sofa::helper::vector<unsigned int> & /*indices*/);

        /// Callback to add tetrahedron elements.
        void ApplyTopologyChange(const core::topology::TetrahedraAdded* /*event*/);
        /// Callback to remove tetrahedron elements.
        void ApplyTopologyChange(const core::topology::TetrahedraRemoved* /*event*/);


        ///////////////////////// Functions on Hexahedron //////////////////////////////////////

        /// Mass coefficient Creation/Destruction functions for Hexahedral Mesh:
        /// Vertex coefficient of mass matrix creation function to handle creation of new hexahedra
        void applyHexahedronCreation(const sofa::helper::vector< unsigned int >& /*indices*/,
                const sofa::helper::vector< Hexahedron >& /*elems*/,
                const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
                const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/);

        /// Vertex coefficient of mass matrix destruction function to handle creation of new hexahedra
        void applyHexahedronDestruction(const sofa::helper::vector<unsigned int> & /*indices*/);

        /// Callback to add hexahedron elements.
        virtual void ApplyTopologyChange(const core::topology::HexahedraAdded* /*event*/);
         /// Callback to remove hexahedron elements.
        virtual void ApplyTopologyChange(const core::topology::HexahedraRemoved* /*event*/);

    protected:
        HighOrderMeshMatrixMass<DataTypes,TMassType>* m;
    };
    VertexMassHandler* vertexMassHandler;

    class EdgeMassHandler : public topology::TopologyDataHandler<Edge,MassVector>
    {
    public:
        EdgeMassHandler(HighOrderMeshMatrixMass<DataTypes,TMassType>* _m, topology::EdgeData<helper::vector<TMassType> >* _data) : topology::TopologyDataHandler<Edge,helper::vector<TMassType> >(_data), m(_m) {}

        /// Edge mass coefficient matrix creation function
        void applyCreateFunction(unsigned int edgeIndex, MassType & EdgeMass,
                const Edge&,
                const sofa::helper::vector< unsigned int > &,
                const sofa::helper::vector< double >&);

        using topology::TopologyDataHandler<Edge,MassVector>::ApplyTopologyChange;

        ///////////////////////// Functions on Triangles //////////////////////////////////////

        /// Edge coefficient of mass matrix creation function to handle creation of new triangles
        void applyTriangleCreation(const sofa::helper::vector< unsigned int >& /*indices*/,
                const sofa::helper::vector< Triangle >& /*elems*/,
                const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
                const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/);

        /// Edge coefficient of mass matrix destruction function to handle creation of new triangles
        void applyTriangleDestruction(const sofa::helper::vector<unsigned int> & /*indices*/);

        /// Callback to add triangles elements.
        void ApplyTopologyChange(const core::topology::TrianglesAdded* /*event*/);
        /// Callback to remove triangles elements.
        void ApplyTopologyChange(const core::topology::TrianglesRemoved* /*event*/);


        ///////////////////////// Functions on Quads //////////////////////////////////////

        /// Edge coefficient of mass matrix creation function to handle creation of new quads
        void applyQuadCreation(const sofa::helper::vector< unsigned int >& /*indices*/,
                const sofa::helper::vector< Quad >& /*elems*/,
                const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
                const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/);

        /// Edge coefficient of mass matrix destruction function to handle creation of new quads
        void applyQuadDestruction(const sofa::helper::vector<unsigned int> & /*indices*/);

        /// Callback to add quads elements.
        void ApplyTopologyChange(const core::topology::QuadsAdded* /*event*/);
        /// Callback to remove quads elements.
        void ApplyTopologyChange(const core::topology::QuadsRemoved* /*event*/);


        ///////////////////////// Functions on Tetrahedron //////////////////////////////////////

        /// Edge coefficient of mass matrix creation function to handle creation of new tetrahedra
        void applyTetrahedronCreation(const sofa::helper::vector< unsigned int >& /*indices*/,
                const sofa::helper::vector< Tetrahedron >& /*elems*/,
                const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
                const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/);

        /// Edge coefficient of mass matrix destruction function to handle creation of new tetrahedra
        void applyTetrahedronDestruction(const sofa::helper::vector<unsigned int> & /*indices*/);

        /// Callback to add tetrahedron elements.
        void ApplyTopologyChange(const core::topology::TetrahedraAdded* /*event*/);
        /// Callback to remove tetrahedron elements.
        void ApplyTopologyChange(const core::topology::TetrahedraRemoved* /*event*/);


        ///////////////////////// Functions on Hexahedron //////////////////////////////////////

        /// Edge coefficient of mass matrix creation function to handle creation of new hexahedra
        void applyHexahedronCreation(const sofa::helper::vector< unsigned int >& /*indices*/,
                const sofa::helper::vector< Hexahedron >& /*elems*/,
                const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
                const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/);

        /// Edge coefficient of mass matrix destruction function to handle creation of new hexahedra
        void applyHexahedronDestruction(const sofa::helper::vector<unsigned int> & /*indices*/);

        /// Callback to add hexahedron elements.
        void ApplyTopologyChange(const core::topology::HexahedraAdded* /*event*/);
         /// Callback to remove hexahedron elements.
        void ApplyTopologyChange(const core::topology::HexahedraRemoved* /*event*/);

    protected:
        HighOrderMeshMatrixMass<DataTypes,TMassType>* m;
    };

    EdgeMassHandler* edgeMassHandler;

    class TetrahedronMassHandler : public topology::TopologyDataHandler<Tetrahedron,MassVectorVector>
    {
    public:
        typedef typename DataTypes::Real Real;
        TetrahedronMassHandler(HighOrderMeshMatrixMass<DataTypes,TMassType>* _m, topology::TetrahedronData<helper::vector<MassVector> >* _data) : topology::TopologyDataHandler<Tetrahedron,helper::vector<MassVector> >(_data), m(_m) {}

        /// Edge mass coefficient matrix creation function
        void applyCreateFunction(unsigned int tetrahedronIndex, MassVector & tetrahedronMass,
                const Tetrahedron&,
                const sofa::helper::vector< unsigned int > &,
                const sofa::helper::vector< double >&);

               /// Edge coefficient of mass matrix destruction function to handle creation of new tetrahedra
//        void applyDestructionFunction(const sofa::helper::vector<unsigned int> & /*indices*/);

    protected:
        HighOrderMeshMatrixMass<DataTypes,TMassType>* m;
    };
	 class TriangleMassHandler : public topology::TopologyDataHandler<Triangle,MassVectorVector>
    {
    public:
        typedef typename DataTypes::Real Real;
        TriangleMassHandler(HighOrderMeshMatrixMass<DataTypes,TMassType>* _m, topology::TriangleData<helper::vector<MassVector> >* _data) : topology::TopologyDataHandler<Triangle,helper::vector<MassVector> >(_data), m(_m) {}

        /// Edge mass coefficient matrix creation function
        void applyCreateFunction(unsigned int triangleIndex, MassVector & triangleMass,
                const Triangle&,
                const sofa::helper::vector< unsigned int > &,
                const sofa::helper::vector< double >&);

               /// Edge coefficient of mass matrix destruction function to handle creation of new tetrahedra
//        void applyDestructionFunction(const sofa::helper::vector<unsigned int> & /*indices*/);

    protected:
        HighOrderMeshMatrixMass<DataTypes,TMassType>* m;
    };

    TetrahedronMassHandler* tetrahedronMassHandler;
    TriangleMassHandler* triangleMassHandler;
	// the array where the weights and coordinates of each integration points are stored
	std::vector<NumericalIntegrationMassData> numericalIntegrationStiffnessDataArray;

	std::vector< std::vector<RegularMassCoefficientEntry> > regularMassCoefficientEntryArray;
	std::vector< std::vector<Real> > bezierRegularMassCoefficientEntryArray;

};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_HIGHORDERDFEM_MESHMATRIXMASS_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_HIGHORDER_FEM_API HighOrderMeshMatrixMass<defaulttype::Vec3dTypes,double>;
extern template class SOFA_HIGHORDER_FEM_API HighOrderMeshMatrixMass<defaulttype::Vec2dTypes,double>;
extern template class SOFA_HIGHORDER_FEM_API HighOrderMeshMatrixMass<defaulttype::Vec1dTypes,double>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_HIGHORDER_FEM_API HighOrderMeshMatrixMass<defaulttype::Vec3fTypes,float>;
extern template class SOFA_HIGHORDER_FEM_API HighOrderMeshMatrixMass<defaulttype::Vec2fTypes,float>;
extern template class SOFA_HIGHORDER_FEM_API HighOrderMeshMatrixMass<defaulttype::Vec1fTypes,float>;
#endif
#endif

} // namespace mass

} // namespace component

} // namespace sofa

#endif
