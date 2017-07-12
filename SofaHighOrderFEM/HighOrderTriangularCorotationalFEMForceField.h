
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



template<class DataTypes>
class SOFA_HIGHORDER_FEM_API HighOrderTriangularCorotationalFEMForceField : public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(HighOrderTriangularCorotationalFEMForceField,DataTypes), SOFA_TEMPLATE(core::behavior::ForceField,DataTypes));

    typedef core::behavior::ForceField<DataTypes> Inherited;
    typedef typename DataTypes::Real        Real        ;
    typedef typename DataTypes::Coord       Coord       ;
    typedef typename DataTypes::Deriv       Deriv       ;
    typedef typename DataTypes::VecCoord    VecCoord    ;
    typedef typename DataTypes::VecDeriv    VecDeriv    ;
    typedef typename DataTypes::VecReal     VecReal     ;
    typedef Data<VecCoord>                  DataVecCoord;
    typedef Data<VecDeriv>                  DataVecDeriv;    
	typedef sofa::core::topology::Topology::Triangle Triangle;
	typedef sofa::core::topology::Topology::TriangleID TriangleID;
	typedef sofa::core::topology::Topology::Point Point;
	typedef sofa::core::topology::Topology::Edge Edge;
	typedef sofa::core::topology::Topology::Quad Quad;
	typedef sofa::core::topology::BaseMeshTopology::EdgesInTriangle EdgesInTriangle;


	typedef helper::vector<Coord> SetAnisotropyDirectionArray; // When the model is anisotropic, for instance in invariant I4
	typedef helper::vector<Real> SetParameterArray; //necessary to store hyperelastic parameters (mu, lambda ...)

    typedef Mat<2,2,Real>       Mat2x2  ;
    typedef Mat<3,3,Real>       Mat3x3  ;
    typedef Mat<3,2,Real>       Mat3x2  ;
    typedef Mat<6,3,Real>       Mat6x3  ;
    typedef Mat<3,4,Real>       Mat3x4  ;

	typedef Mat<15,3,Real>		Mat15x3  ;
	typedef Mat<15,4,Real>		Mat15x4  ;
	typedef Mat<45,3,Real>     Mat45x3  ;
	typedef Mat<45,4,Real>		Mat45x4  ;
	typedef Mat<105,3,Real>     Mat105x3  ;
	typedef Mat<105,4,Real>		Mat105x4  ;
	typedef Mat<210,3,Real>    Mat210x3  ;
	typedef Mat<210,4,Real>		Mat210x4  ;


	typedef typename defaulttype::MatSym<3,Real> MatrixSym;
    // In case of non 2D template
    typedef Vec<2,Real> Vec2;
    typedef Vec<3,Real> Vec3;
    typedef Vec<4,Real> Vec4;

    typedef StdVectorTypes< Vec2, Vec2, Real >     GeometricalTypes ; /// assumes the geometry object type is 2D
	typedef typename sofa::component::topology::HighOrderTriangleSetGeometryAlgorithms<GeometricalTypes>::VecPointID VecPointID;

    typedef helper::vector<Real> ParameterArray;
    typedef helper::vector<Coord> AnisotropyDirectionArray;

    typedef enum
    {
        POLAR_DECOMPOSITION,
        QR_DECOMPOSITION,
	    POLAR_DECOMPOSITION_MODIFIED,
		LINEAR_ELASTIC
    } RotationDecompositionMethod;

	/// the way the stiffness matrix should be computed on HighOrder elements
	typedef enum 
	{
        AFFINE_ELEMENT_INTEGRATION = 1,
        NUMERICAL_INTEGRATION = 2,
        STANDARD_INTEGRATION = 3,
        BEZIER_NUMERICAL_INTEGRATION = 4,
        NUMERICAL_INTEGRATION_2 = 5
	} IntegrationMethod;

		/// the way the stiffness matrix should be computed on HighOrder elements
	typedef enum 
	{
		ISOTROPIC=1,
		TRANSVERSE_ISOTROPIC=2,
		ORTHOTROPIC=3,
		CUBIC=4
	} ElasticitySymmetry;

protected:
		// structure that store coefficients matrices for elements of order 2,3,4,5
	// to save memory only store pointer to matrices of predefined size
	struct weightArrayPointer {
	public:
		boost::shared_ptr<Mat15x3>   weightArrayQuadratic[2];
		boost::shared_ptr<Mat45x3>  weightArrayCubic[2];
		boost::shared_ptr<Mat105x3>  weightArrayQuartic[2];
		boost::shared_ptr<Mat210x3> weightArrayQuintic[2];

		void allocate(size_t degree);
	};
	// the array where stiffness coefficients are stored for affine elements of any degree
	std::vector<Vec3> affineStiffnessCoefficientArray;
	// the array where stiffness coefficients are stored for affine elements of degree < 5 
	weightArrayPointer affineStiffnessCoefficientPreStoredArray;
		// the data stored for each integration point
	struct NumericalIntegrationStiffnessData {
		// the weight of the integration point
		Real integrationWeight;
		// for each pair of control point store the weight matrix 6 w_\gamma * dN_p/d\param_i(\param_gamma) *  dN_q/d\param_j (\param_gamma)
		std::vector<Vec3> weightArray;
		// for each pair of control point store the weight matrix 6 w_\gamma * dN_p/d\param_i(\param_gamma) *  dN_q/d\param_j (\param_gamma)
		weightArrayPointer arrayPointer;
		// for each control point  store the derivative of the shape functions
		std::vector<Deriv> coefficientArray;
		/// barycentric coordinate of the integration point \param_i
		Vec3 integrationPoint;
	};
	// the array where the weights and coordinates of each integration points are stored
	std::vector<NumericalIntegrationStiffnessData> numericalIntegrationStiffnessDataArray;
	/// the elasticity tensor 
	Mat3x3 elasticityTensor;

    /// data structure stored for each triangle
    class TriangleRestInformation
	{
	public:
		typedef typename DataTypes::Real  Real;
		typedef Mat<2,2,Real> Mat3x3;
		

		/// rest volume
      
		Coord shapeVector[3];
        Coord restEdgeVector[3];
        Mat3x3 rotation; // rotation from deformed to rest configuration
        Mat3x3 restRotation; // used for QR decomposition
		helper::vector<Mat3x3> stiffnessVector; // the nc*(nc+1)/2 stiffness matrices where nc is the number of control points
		size_t v[3]; // the indices of the 4 vertices
		helper::vector<Mat3x3> rotatedStiffnessVector; // the nc*(nc+1)/2 stiffness matrices where nc is the number of control points
		helper::vector<Mat3x3> reducedStiffnessVector;
		// store 6 2x3 matrices per integration points
		helper::vector<  Mat3x4 >  integrationPointsStiffnessVector;
		// store the 3 rest edge vector for each integration point
		helper::vector< helper::vector<Coord> > integrationPointsRestEdgeVector;
		// store a rest rotation for each integration point
		helper::vector< Mat2x2 > integrationPointsRestRotationArray;

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
        typedef typename HighOrderTriangularCorotationalFEMForceField<DataTypes>::TriangleRestInformation TriangleRestInformation;

        FTCFTriangleHandler(HighOrderTriangularCorotationalFEMForceField<DataTypes>* ff,
                TriangleData<sofa::helper::vector<TriangleRestInformation> >* data )
            :TopologyDataHandler<Triangle, sofa::helper::vector<TriangleRestInformation> >(data)
            ,ff(ff)
        {

        }

        void applyCreateFunction(unsigned int, TriangleRestInformation &t, const Triangle
                &, const sofa::helper::vector<unsigned int> &, const sofa::helper::vector<double> &);

    protected:
        HighOrderTriangularCorotationalFEMForceField<DataTypes>* ff;

    };

    TriangleData<sofa::helper::vector<TriangleRestInformation> > triangleInfo;


    sofa::core::topology::BaseMeshTopology* _topology;
	sofa::component::topology::HighOrderTriangleSetGeometryAlgorithms<GeometricalTypes>* highOrderTriangleGeo;
    VecCoord  _initialPoints;///< the intial positions of the points

    bool updateMatrix;
    bool updateTopologyInfo;


    RotationDecompositionMethod decompositionMethod;

	ElasticitySymmetry elasticitySymmetry;

	std::vector<Mat2x2> anisotropyMatrixArray;
	std::vector<Real> anisotropyScalarArray;

    Real lambda;  /// first Lame coefficient
    Real mu;    /// second Lame coefficient

	
	IntegrationMethod    integrationMethod;



    HighOrderTriangularCorotationalFEMForceField();

    virtual ~HighOrderTriangularCorotationalFEMForceField();

public:
	// user input
	Data<std::string> d_method; ///< the computation method of the displacements


	Data<Real> d_poissonRatio; // stiffness coefficient for isotropic elasticity;
	Data<Real> d_youngModulus;

	Data<std::string> d_anisotropy; // the type of isotropy
	Data<ParameterArray> d_anisotropyParameter; // the set of parameters defining the elasticity anisotropy
	Data<AnisotropyDirectionArray> d_anisotropyDirection; // the directions of anisotropy

	/// the order of integration for numerical integration
	Data<size_t>	     numericalIntegrationOrder;
	/// the type of numerical integration method chosen
	Data<std::string>	     numericalIntegrationMethod;
	/// the type of integration method chosen for non linear element.
	Data<std::string>	 d_integrationMethod; 
	// if one rotation is attached to each integration point 
	Data<bool> d_oneRotationPerIntegrationPoint;
	// measure the time spent in assembling the stiffness matrix
	Data<Real> d_assemblyTime;
	// whether each affine element should be assembled with the affine assembly method irrespective to the chosen integration method  
	Data<bool> d_forceAffineAssemblyForAffineElements;

    virtual void init();


    virtual void addForce(const sofa::core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv &  dataF, const DataVecCoord &  dataX , const DataVecDeriv & dataV ) ;
    virtual void addDForce(const sofa::core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, DataVecDeriv&   datadF , const DataVecDeriv&   datadX ) ;
    virtual SReal getPotentialEnergy(const core::MechanicalParams*, const DataVecCoord&) const;

    void updateTopologyInformation();

    virtual Real getLambda() const { return lambda;}
    virtual Real getMu() const { return mu;}

    void setYoungModulus(const double modulus)
    {
        d_youngModulus.setValue((Real)modulus);
    }
    void setPoissonRatio(const double ratio)
    {
        d_poissonRatio.setValue((Real)ratio);
    }
    void setRotationDecompositionMethod( const RotationDecompositionMethod m)
    {
        decompositionMethod=m;
    }
    void draw(const core::visual::VisualParams* vparams);
    /// compute lambda and mu based on the Young modulus and Poisson ratio
    void updateLameCoefficients();


	void computeTriangleStiffnessEdgeMatrix(const Coord position[3],Mat3x4 edgeStiffness[2]);
	friend class FTCFTriangleHandler;

protected :
    FTCFTriangleHandler* triangleHandler;
	// for profiling the assembly task
	helper::system::thread::ctime_t totalUpdateMat;
	static void computeQRRotation( Mat2x2 &r, const Coord *dp);

	virtual const helper::vector<Mat2x2> &getStiffnessArray(const size_t i,
		TriangleRestInformation *restTriangle);

	virtual const helper::vector<Mat2x2> &getRotatedStiffnessArray(const size_t i,
		const TriangleRestInformation *restTriangle);
	// update anisotropyMatrixArray and anisotropyScalarArray based on input data
	void assembleAnisotropicTensors();
	// compute elasticity tensor for isotropic and anisotropic cases
	void computeElasticityTensor(); 		

};


} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_FASTTRIANGULARCOROTATIONALFORCEFIELD_H
