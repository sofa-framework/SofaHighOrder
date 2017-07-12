
#ifndef SOFA_COMPONENT_FORCEFIELD_HIGHORDERTETRAHEDRALCOROTATIONALFEMFORCEFIELD_H
#define SOFA_COMPONENT_FORCEFIELD_HIGHORDERTETRAHEDRALCOROTATIONALFEMFORCEFIELD_H

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif
#include "initHighOrderFEM.h"
#include <sofa/core/behavior/ForceField.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/helper/fixed_array.h>
#include <sofa/helper/vector.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>
#include <SofaBaseTopology/TopologyData.h>
#include <sofa/defaulttype/MatSym.h>
#include <sofa/helper/system/thread/CTime.h>

namespace sofa
{

namespace component
{
namespace topology
{
	template< class DataTypes> class HighOrderTetrahedronSetGeometryAlgorithms;
}
namespace forcefield
{

using namespace sofa::defaulttype;
using namespace sofa::component::topology;
using namespace sofa::component::container;


template<class DataTypes>
class SOFA_HIGHORDER_FEM_API HighOrderTetrahedralDiffusionForceField : public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(HighOrderTetrahedralDiffusionForceField,DataTypes), SOFA_TEMPLATE(core::behavior::ForceField,DataTypes));

    typedef core::behavior::ForceField<DataTypes> Inherited;
    typedef typename DataTypes::Real        Real        ;
    typedef typename DataTypes::Coord       Coord       ;
    typedef typename DataTypes::Deriv       Deriv       ;
    typedef typename DataTypes::VecCoord    VecCoord    ;
    typedef typename DataTypes::VecDeriv    VecDeriv    ;
    typedef Data<VecCoord>                  DataVecCoord;
    typedef Data<VecDeriv>                  DataVecDeriv;    
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
	typedef sofa::core::topology::BaseMeshTopology::TetrahedraAroundTriangle TetrahedraAroundTriangle;

	typedef helper::vector<Coord> SetAnisotropyDirectionArray; // When the model is anisotropic, for instance in invariant I4
	typedef helper::vector<Real> SetParameterArray; //necessary to store hyperelastic parameters (mu, lambda ...)

	typedef Mat<3,3,Real>       Mat3x3  ;
	typedef Mat<4,4,Real>       Mat4x4  ;
	typedef Mat<45,6,Real>		Mat45x6  ;
	typedef Vec<45,Real>		Vec45  ;
	typedef Mat<190,6,Real>     Mat190x6  ;
	typedef Vec<190,Real>		Vec190  ;
	typedef Mat<595,6,Real>     Mat595x6  ;
	typedef Vec<595,Real>		Vec595  ;
	typedef Mat<1540,6,Real>    Mat1540x6  ;
	typedef Vec<1540,Real>		Vec1540  ;

    // In case of non 3D template
    typedef Vec<3,Real> Vec3;
    typedef Vec<4,Real> Vec4;
	typedef Vec<6,Real> Vec6;


	// Vectors GPU compatible
	typedef typename VecCoord::template rebind<unsigned int>::other VecUInt;
	typedef typename VecCoord::template rebind<Real>::other VecReal;

	/// assumes the mechanical object type (3D)
	typedef StdVectorTypes< Vec3, Vec3, Real >     MechanicalTypes ;
	typedef MechanicalObject<MechanicalTypes>      MechObject;
	typedef typename sofa::component::topology::HighOrderTetrahedronSetGeometryAlgorithms<MechanicalTypes>::VecPointID VecPointID;

    typedef helper::vector<Real> ParameterArray;
    typedef helper::vector<size_t> DOFArray;
    typedef helper::vector<Vec3> AnisotropyDirectionArray;


	/// the way the stiffness matrix should be computed on HighOrder elements
	typedef enum 
	{
		AFFINE_ELEMENT_INTEGRATION=1,
		NUMERICAL_INTEGRATION=2,
		STANDARD_INTEGRATION=3,
        BEZIER_NUMERICAL_INTEGRATION = 4,
		NUMERICAL_INTEGRATION_2=5
	} IntegrationMethod;

		/// the way the stiffness matrix should be computed on HighOrder elements
	typedef enum 
	{
		ISOTROPIC=1,
		TRANSVERSE_ISOTROPIC=2,
		ORTHOTROPIC=3
	} DiffusionSymmetry;

protected:
	// structure that store coefficients matrices for elements of order 2,3,4,5
	// to save memory only store pointer to matrices of predefined size
	struct weightArrayPointer {
	public:
		boost::shared_ptr<Mat45x6>   weightArrayQuadratic;
		boost::shared_ptr<Mat190x6>  weightArrayCubic;
		boost::shared_ptr<Mat595x6>  weightArrayQuartic;
		boost::shared_ptr<Mat1540x6> weightArrayQuintic;

		boost::shared_ptr<Vec45>   resultQuadratic;
		boost::shared_ptr<Vec190>  resultCubic;
		boost::shared_ptr<Vec595>  resultQuartic;
		boost::shared_ptr<Vec1540> resultQuintic;

		void allocate(size_t degree);
	};
	struct weightVectorPointer {
	public:

		boost::shared_ptr<Vec45>   resultQuadratic;
		boost::shared_ptr<Vec190>  resultCubic;
		boost::shared_ptr<Vec595>  resultQuartic;
		boost::shared_ptr<Vec1540> resultQuintic;

		void allocate(size_t degree);
		void resetResult(size_t degree);
		void updateResult(size_t degree,weightArrayPointer & wap,Vec6 &input);
		void updateArray(size_t degree,sofa::helper::vector<Real> &array);
	};
	// the array where stiffness coefficients are stored for affine elements.
	std::vector<Vec6> affineStiffnessCoefficientArray;
	// the array where stiffness coefficients are stored for affine elements of degree < 5 
	weightArrayPointer affineStiffnessCoefficientPreStoredArray;
	// the data stored for each integration point
	struct NumericalIntegrationStiffnessData {
		// the weight of the integration point
		Real integrationWeight;
		// for each pair of control point store the weight matrix 6 w_\gamma * dN_p/d\param_i(\param_gamma) *  dN_q/d\param_j (\param_gamma)
		std::vector<Mat4x4> weightArray;
		// for each pair of control point store the weight matrix 6 w_\gamma * dN_p/d\param_i(\param_gamma) *  dN_q/d\param_j (\param_gamma)
		std::vector<Vec6> weightVectorizedArray;
		// for each control point  store the derivative of the shape functions
		std::vector<Vec3> coefficientArray;
		weightArrayPointer arrayPointer;
		/// barycentric coordinate of the integration point \param_i
		Vec4 integrationPoint;
	};
	// the array where the weights and coordinates of each integration points are stored
	std::vector<NumericalIntegrationStiffnessData> numericalIntegrationStiffnessDataArray;
	/// the elasticity tensor 
	Mat3x3 diffusionTensor;

    /// data structure stored for each tetrahedron
    class TetrahedronRestInformation
	{
	public:
		typedef typename DataTypes::Real  Real;
		typedef Mat<3,3,Real> Mat3x3;
		

		/// rest volume
     
		helper::vector<Real> stiffnessVector; // the nc*(nc+1)/2 stiffness matrices where nc is the number of control points
		size_t v[4]; // the indices of the 4 vertices

		// store 6 2x3 matrices per integration points
		helper::vector<  helper::vector<Real> >  integrationPointsStiffnessVector;

		// store a rest rotation for each integration point
		helper::vector< Mat3x3 > integrationPointsRestRotationArray;

        /// Output stream
        inline friend std::ostream& operator<< ( std::ostream& os, const TetrahedronRestInformation& /*eri*/ )
        {
            return os;
        }

        /// Input stream
        inline friend std::istream& operator>> ( std::istream& in, TetrahedronRestInformation& /*eri*/ )
        {
            return in;
        }

        TetrahedronRestInformation()
        {
        }
    };

    class FTCFTetrahedronHandler : public TopologyDataHandler<Tetrahedron, sofa::helper::vector<TetrahedronRestInformation> >
    {
    public:
        typedef typename HighOrderTetrahedralDiffusionForceField<DataTypes>::TetrahedronRestInformation TetrahedronRestInformation;

        FTCFTetrahedronHandler(HighOrderTetrahedralDiffusionForceField<DataTypes>* ff,
                TetrahedronData<sofa::helper::vector<TetrahedronRestInformation> >* data )
            :TopologyDataHandler<Tetrahedron, sofa::helper::vector<TetrahedronRestInformation> >(data)
            ,ff(ff)
        {

        }

        void applyCreateFunction(unsigned int, TetrahedronRestInformation &t, const Tetrahedron
                &, const sofa::helper::vector<unsigned int> &, const sofa::helper::vector<double> &);

    protected:
        HighOrderTetrahedralDiffusionForceField<DataTypes>* ff;

    };

    TetrahedronData<sofa::helper::vector<TetrahedronRestInformation> > tetrahedronInfo;

    sofa::core::topology::BaseMeshTopology* _topology;
	/// Pointer to mechanical mechanicalObject
	MechanicalObject<MechanicalTypes> *mechanicalObject;
	sofa::component::topology::HighOrderTetrahedronSetGeometryAlgorithms<MechanicalTypes>* highOrderTetraGeo;


    bool updateMatrix;
    bool updateTopologyInfo;
	 
	DiffusionSymmetry diffusionSymmetry;

	IntegrationMethod    integrationMethod;

    HighOrderTetrahedralDiffusionForceField();

    virtual ~HighOrderTetrahedralDiffusionForceField();

public:
	// data for user input

	Data<std::string> d_anisotropy; // the type of isotropy
	/// the order of integration for numerical integration
	Data<size_t>	     numericalIntegrationOrder;
	/// the type of numerical integration method chosen
	Data<std::string>	     numericalIntegrationMethod;
	Data<Real> d_diffusivity; // stiffness coefficient for isotropic elasticity;
	/// anisotropy parameters and directions
	Data<ParameterArray> d_anisotropyParameter;
	Data<AnisotropyDirectionArray> d_anisotropyDirection;
	/// specify the set of indices of the state vector where diffusion should be performed
	Data<DOFArray> m_diffusionDOF;
	/// Mechanic xml tags of the system.
	Data<std::string> m_tagMeshMechanics;
	/// the type of integration method chosen for non linear element.
	Data<std::string>	 d_integrationMethod; 
	// measure the time spent in assembling the stiffness matrix
	Data<Real> d_assemblyTime;
	// whether each affine element should be assembled with the affine assembly method irrespective to the chosen integration method  
	Data<bool> d_forceAffineAssemblyForAffineElements;
	virtual void init();


    virtual void addForce(const sofa::core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv &  dataF, const DataVecCoord &  dataX , const DataVecDeriv & dataV ) ;
    virtual void addDForce(const sofa::core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv&   datadF , const DataVecDeriv&   datadX ) ;
    virtual SReal getPotentialEnergy(const core::MechanicalParams*, const DataVecCoord&) const;

    void updateTopologyInformation();


    void draw(const core::visual::VisualParams* vparams);
    /// compute lambda and mu based on the Young modulus and Poisson ratio
    void updateLameCoefficients();

	void computeTetrahedronStiffnessEdgeMatrix(const Vec3 position[4],Vec6 &edgeStiffness);

	friend class FTCFTetrahedronHandler;

protected :
    FTCFTetrahedronHandler* tetrahedronHandler;
	// for profiling the assembly task
	helper::system::thread::ctime_t totalUpdateMat;

	virtual const helper::vector<Real> &getStiffnessArray(const size_t i,
		const TetrahedronRestInformation *restTetra);

	// compute elasticity tensor for isotropic and anisotropic cases
	void computeDiffusivityTensor(); 		

};


} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_FASTTETRAHEDRALCOROTATIONALFORCEFIELD_H
