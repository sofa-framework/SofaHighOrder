#ifndef SOFA_COMPONENT_FORCEFIELD_HIGHORDERTRIANGULARCOROTATIONALFEMFORCEFIELD_H
#define SOFA_COMPONENT_FORCEFIELD_HIGHORDERTRIANGULARCOROTATIONALFEMFORCEFIELD_H

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
	template< class DataTypes> class HighOrderTriangleSetGeometryAlgorithms;
}
namespace forcefield
{

using namespace sofa::defaulttype;
using namespace sofa::component::topology;
using namespace sofa::component::container;


template<class DataTypes>
class SOFA_HIGHORDER_FEM_API HighOrderTriangularDiffusionForceField : public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(HighOrderTriangularDiffusionForceField,DataTypes), SOFA_TEMPLATE(core::behavior::ForceField,DataTypes));

    typedef core::behavior::ForceField<DataTypes> Inherited;
    typedef typename DataTypes::Real        Real        ;
    typedef typename DataTypes::Coord       Coord       ;
    typedef typename DataTypes::Deriv       Deriv       ;
    typedef typename DataTypes::VecCoord    VecCoord    ;
    typedef typename DataTypes::VecDeriv    VecDeriv    ;
    typedef Data<VecCoord>                  DataVecCoord;
    typedef Data<VecDeriv>                  DataVecDeriv;    
	typedef sofa::core::topology::Topology::Triangle Triangle;
	typedef sofa::core::topology::Topology::TetraID TetraID;
	typedef sofa::core::topology::Topology::Tetra Tetra;
	typedef sofa::core::topology::Topology::Point Point;
	typedef sofa::core::topology::Topology::Edge Edge;
	typedef sofa::core::topology::Topology::Quad Quad;
	typedef sofa::core::topology::Topology::Hexahedron Hexahedron;
	typedef sofa::core::topology::BaseMeshTopology::EdgesInTriangle EdgesInTriangle;
	typedef sofa::core::topology::BaseMeshTopology::EdgesInQuad EdgesInQuad;
	typedef sofa::core::topology::BaseMeshTopology::EdgesInHexahedron EdgesInHexahedron;
	typedef sofa::core::topology::BaseMeshTopology::TetrahedraAroundTriangle TetrahedraAroundTriangle;

	typedef helper::vector<Coord> SetAnisotropyDirectionArray; // When the model is anisotropic, for instance in invariant I4
	typedef helper::vector<Real> SetParameterArray; //necessary to store hyperelastic parameters (mu, lambda ...)
	typedef Mat<2,2,Real>       Mat2x2  ;
	typedef Mat<3,3,Real>       Mat3x3  ;

	typedef Mat<15,3,Real>		Mat15x3  ;
	typedef Vec<15,Real>		Vec15  ;
	typedef Mat<45,3,Real>     Mat45x3  ;
	typedef Vec<45,Real>		Vec45  ;
	typedef Mat<105,3,Real>     Mat105x3  ;
	typedef Vec<105,Real>		Vec105  ;
	typedef Mat<210,3,Real>    Mat210x3  ;
	typedef Vec<210,Real>		Vec210  ;

    // In case of non 3D template
	typedef Vec<2,Real> Vec2;
	typedef Vec<3,Real> Vec3;



	// Vectors GPU compatible
	typedef typename VecCoord::template rebind<unsigned int>::other VecUInt;
	typedef typename VecCoord::template rebind<Real>::other VecReal;

	/// assumes the mechanical object type (3D)
	typedef StdVectorTypes< Vec2, Vec2, Real >     MechanicalTypes ;
	typedef MechanicalObject<MechanicalTypes>      MechObject;
	typedef typename sofa::component::topology::HighOrderTriangleSetGeometryAlgorithms<MechanicalTypes>::VecPointID VecPointID;

    typedef helper::vector<Real> ParameterArray;
    typedef helper::vector<size_t> DOFArray;
    typedef helper::vector<Vec2> AnisotropyDirectionArray;


	/// the way the stiffness matrix should be computed on HighOrder elements
	typedef enum 
	{
		AFFINE_ELEMENT_INTEGRATION=1,
		NUMERICAL_INTEGRATION=2,
		STANDARD_INTEGRATION=3,
		NUMERICAL_INTEGRATION_2=4
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
		boost::shared_ptr<Mat15x3>   weightArrayQuadratic;
		boost::shared_ptr<Mat45x3>  weightArrayCubic;
		boost::shared_ptr<Mat105x3>  weightArrayQuartic;
		boost::shared_ptr<Mat210x3> weightArrayQuintic;

		void allocate(size_t degree);
	};
	// the array where stiffness coefficients are stored for affine elements.
	std::vector<Vec3> affineStiffnessCoefficientArray;
	// the array where stiffness coefficients are stored for affine elements of degree < 5 
	weightArrayPointer affineStiffnessCoefficientPreStoredArray;
	// the data stored for each integration point
	struct NumericalIntegrationStiffnessData {
		// the weight of the integration point
		Real integrationWeight;
		// for each pair of control point store the weight matrix 6 w_\gamma * dN_p/d\param_i(\param_gamma) *  dN_q/d\param_j (\param_gamma)
		std::vector<Mat3x3> weightArray;
		// for each pair of control point store the weight matrix 6 w_\gamma * dN_p/d\param_i(\param_gamma) *  dN_q/d\param_j (\param_gamma)
		std::vector<Vec3> weightVectorizedArray;
		// for each control point  store the derivative of the shape functions
		std::vector<Vec2> coefficientArray;
		weightArrayPointer arrayPointer;
		/// barycentric coordinate of the integration point \param_i
		Vec3 integrationPoint;
	};
	// the array where the weights and coordinates of each integration points are stored
	std::vector<NumericalIntegrationStiffnessData> numericalIntegrationStiffnessDataArray;
	/// the elasticity tensor 
	Mat2x2 diffusionTensor;

    /// data structure stored for each triangle
    class TriangleRestInformation
	{
	public:
		typedef typename DataTypes::Real  Real;
		typedef Mat<3,3,Real> Mat3x3;
		

		/// rest volume
     
		helper::vector<Real> stiffnessVector; // the nc*(nc+1)/2 stiffness matrices where nc is the number of control points
		size_t v[3]; // the indices of the 3 vertices

        /// Output stream
        inline friend std::ostream& operator<< ( std::ostream& os, const TriangleRestInformation& /*eri*/ )
        {
            return os;
        }

        /// Input stream
        inline friend std::istream& operator>> ( std::istream& in, TriangleRestInformation& /*eri*/ )
        {
            return in;
        }

        TriangleRestInformation()
        {
        }
    };

    class FTCFTriangleHandler : public TopologyDataHandler<Triangle, sofa::helper::vector<TriangleRestInformation> >
    {
    public:
        typedef typename HighOrderTriangularDiffusionForceField<DataTypes>::TriangleRestInformation TriangleRestInformation;

        FTCFTriangleHandler(HighOrderTriangularDiffusionForceField<DataTypes>* ff,
                TriangleData<sofa::helper::vector<TriangleRestInformation> >* data )
            :TopologyDataHandler<Triangle, sofa::helper::vector<TriangleRestInformation> >(data)
            ,ff(ff)
        {

        }

        void applyCreateFunction(unsigned int, TriangleRestInformation &t, const Triangle
                &, const sofa::helper::vector<unsigned int> &, const sofa::helper::vector<double> &);

    protected:
        HighOrderTriangularDiffusionForceField<DataTypes>* ff;

    };

    TriangleData<sofa::helper::vector<TriangleRestInformation> > triangleInfo;

    sofa::core::topology::BaseMeshTopology* _topology;
	/// Pointer to mechanical mechanicalObject
	MechanicalObject<MechanicalTypes> *mechanicalObject;
	sofa::component::topology::HighOrderTriangleSetGeometryAlgorithms<MechanicalTypes>* highOrderTrianGeo;


    bool updateMatrix;
    bool updateTopologyInfo;
	 
	DiffusionSymmetry diffusionSymmetry;

	IntegrationMethod    integrationMethod;

    HighOrderTriangularDiffusionForceField();

    virtual ~HighOrderTriangularDiffusionForceField();

public:
	// data for user input

	Data<std::string> d_anisotropy; // the type of isotropy
	/// the order of integration for numerical integration
	Data<size_t>	     numericalIntegrationOrder;
	/// the type of numerical integration method chosen
	Data<size_t>	     numericalIntegrationMethod;
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

	void computeTriangleStiffnessEdgeMatrix(const Vec2 position[3],Vec3 &edgeStiffness);

	friend class FTCFTriangleHandler;

protected :
    FTCFTriangleHandler* triangleHandler;
	// for profiling the assembly task
	helper::system::thread::ctime_t totalUpdateMat;

	virtual const helper::vector<Real> &getStiffnessArray(const size_t i,
		const TriangleRestInformation *restTrian);

	// compute elasticity tensor for isotropic and anisotropic cases
	void computeDiffusivityTensor(); 		

};


} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_FASTTRIANGULARCOROTATIONALFORCEFIELD_H
