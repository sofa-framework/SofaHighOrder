#include <SofaTest/Sofa_test.h>
#include<sofa/helper/system/SetDirectory.h>
#include <sofa/helper/system/FileRepository.h>
#include <sofa/core/ExecParams.h>

//Including Simulation
#include <sofa/simulation/Simulation.h>
#include <SofaSimulationGraph/DAGSimulation.h>
#include <sofa/simulation/Node.h>

// Including constraint, force and mass
#include <../SofaHighOrderTopology/HighOrderTetrahedronSetTopologyContainer.h>
#include <../SofaHighOrderTopology/BezierTetrahedronSetGeometryAlgorithms.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/core/MechanicalParams.h>
#include <SceneCreator/SceneCreator.h>
//#include <../SofaHighOrderFEM//MeshMatrixMass.h>
#include <../SofaHighOrderTopology/GenerateBezierCylinder.h>
#include <../SofaHighOrderTopology/Mesh2HighOrderTopologicalMapping.h>
#include <sofa/defaulttype/VecTypes.h>

namespace sofa {

using std::cout;
using std::cerr;
using std::endl;
using namespace component;
using namespace defaulttype;

/**  Test the topology and geometry of Bezier Tetrahedral meshes. 
No deformation is applied but only the init() function is called to create a Bezier Tetrahedral mesh from a a regular tetraheral mesh */

template <typename _DataTypes>
struct BezierTetrahedronTopology_test : public Sofa_test<typename _DataTypes::Real>
{
    typedef _DataTypes DataTypes;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Deriv Deriv;
	typedef typename DataTypes::Coord Coord;
	typedef typename DataTypes::Real Real;
    typedef typename sofa::component::topology::HighOrderTetrahedronSetTopologyContainer HighOrderTetrahedronSetTopologyContainer;
	typedef typename sofa::component::topology::BezierTetrahedronSetGeometryAlgorithms<DataTypes> BezierTetrahedronSetGeometryAlgorithms;
	typedef typename sofa::component::topology::HighOrderDegreeType BezierDegreeType;
	typedef typename sofa::component::topology::TetrahedronIndexVector TetrahedronBezierIndex;
	typedef typename HighOrderTetrahedronSetTopologyContainer::HighOrderTetrahedronPointLocation BezierTetrahedronPointLocation;
    typedef typename container::MechanicalObject<DataTypes> MechanicalObject;
//    typedef typename sofa::component::mass::MeshMatrixMass<DataTypes,Real>  MeshMatrixMass;

    
    /// Root of the scene graph
    simulation::Node::SPtr root;      
    /// Simulation
    simulation::Simulation* simulation;  

    // Create the context for the scene
    void SetUp()
    { 
        // Init simulation
        sofa::simulation::setSimulation(simulation = new sofa::simulation::graph::DAGSimulation());

         root = simulation::getSimulation()->createNewGraph("root");
    }
	 // Load the scene BezierTetrahedronTopology.scn from the scenes directory
    void loadScene(std::string sceneName)
    {
        // Load the scene from the xml file
	std::string fileName = std::string(SOFAHIGHORDERTOPOLOGY_TEST_SCENES_DIR) + "/" + sceneName;
        root = sofa::simulation::getSimulation()->load(fileName.c_str());
    }
	 // Load the scene BezierTetrahedronTopology.sc from the Scenes directory
    void createScene()
    {
		// GenerateBezierCylinder object
        typename sofa::component::engine::GenerateBezierCylinder<DataTypes>::SPtr eng= sofa::modeling::addNew<sofa::component::engine::GenerateBezierCylinder<DataTypes> >(root,"cylinder");
		eng->f_radius=0.2;
		eng->f_height=1.0;
		// TetrahedronSetTopologyContainer object
		sofa::component::topology::TetrahedronSetTopologyContainer::SPtr container1= sofa::modeling::addNew<sofa::component::topology::TetrahedronSetTopologyContainer>(root,"Container1");
		sofa::modeling::setDataLink(&eng->f_tetrahedra,&container1->d_tetrahedron);
		sofa::modeling::setDataLink(&eng->f_outputTetrahedraPositions,&container1->d_initPoints);
		// TetrahedronSetGeometryAlgorithms object
        typename sofa::component::topology::TetrahedronSetGeometryAlgorithms<DataTypes>::SPtr geo1= sofa::modeling::addNew<sofa::component::topology::TetrahedronSetGeometryAlgorithms<DataTypes> >(root);
		// mechanicalObject object
        typename MechanicalObject::SPtr meca1= sofa::modeling::addNew<MechanicalObject>(root);
		sofa::modeling::setDataLink(&eng->f_outputTetrahedraPositions,&meca1->x);
		// subnode
	    simulation::Node::SPtr bezierNode = root->createChild("BezierTetrahedronTopology");
		// HighOrderTetrahedronSetTopologyContainer
		sofa::component::topology::HighOrderTetrahedronSetTopologyContainer::SPtr container2= sofa::modeling::addNew<sofa::component::topology::HighOrderTetrahedronSetTopologyContainer>(bezierNode,"Container2");
		// Mesh2HighOrderTopologicalMapping
		sofa::component::topology::Mesh2HighOrderTopologicalMapping::SPtr mapping= sofa::modeling::addNew<sofa::component::topology::Mesh2HighOrderTopologicalMapping>(bezierNode,"Mapping");
		mapping->setTopologies(container1.get(),container2.get());
		mapping->bezierTetrahedronDegree=3;
		// mechanicalObject object
        typename MechanicalObject::SPtr meca2= sofa::modeling::addNew<MechanicalObject>(bezierNode,"BezierMechanicalObject");
		// BezierTetrahedronSetGeometryAlgorithms
        typename sofa::component::topology::BezierTetrahedronSetGeometryAlgorithms<DataTypes>::SPtr geo2= sofa::modeling::addNew<sofa::component::topology::BezierTetrahedronSetGeometryAlgorithms<DataTypes> >(bezierNode);
		// MeshMatrixMass
  //      typename MeshMatrixMass::SPtr mass= sofa::modeling::addNew<MeshMatrixMass >(bezierNode,"BezierMass");
//		mass->m_massDensity=1.0;
//		mass->d_integrationMethod.setValue(std::string("analytical"));
	}
	bool testBezierTetrahedronTopology()
	{
		// Init simulation
		sofa::simulation::getSimulation()->init(root.get());
		HighOrderTetrahedronSetTopologyContainer *container=root->get<HighOrderTetrahedronSetTopologyContainer>(root->SearchDown);
		size_t nTetras,elem;
		BezierDegreeType degree=container->getDegree();
		// check the total number of vertices.
		size_t nbPoints=container->getNumberOfTetrahedralPoints()+container->getNumberOfEdges()*(degree-1)+container->getNumberOfTriangles()*(degree-1)*(degree-2)/2+container->getNumberOfTetrahedra()*((degree-1)*(degree-2)*(degree-3)/6);
        if((size_t)container->getNbPoints()!=nbPoints) {
			ADD_FAILURE() << "wrong number of points " <<container->getNbPoints() << " is wrong. It should be  " <<nbPoints  << std::endl;
			return false;
		}

		sofa::helper::vector<TetrahedronBezierIndex> tbiArray=container->getTetrahedronIndexArray();
		
		BezierTetrahedronPointLocation location; 
		size_t elementIndex, elementOffset/*,localIndex*/;
		for (nTetras=0;nTetras<container->getNumberOfTetrahedra();++nTetras) {
			
  			const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=container->getGlobalIndexArrayOfControlPoints(nTetras);
			// check the number of control points per tetrahedron is correct
			nbPoints=(4+6*(degree-1)+2*(degree-1)*(degree-2)+(degree-1)*(degree-2)*(degree-3)/6);
			if(indexArray.size()!=nbPoints) {
				ADD_FAILURE() << "wrong number of control points in tetrahedron " <<nTetras<< ". It is "<<indexArray.size() <<" and should be "<<nbPoints  << std::endl;
				return false;
			}
			for(elem=0;elem<indexArray.size();++elem) {
   				size_t globalIndex=container->getGlobalIndexOfControlPoint(nTetras,tbiArray[elem]);
				// check that getGlobalIndexOfBezierPoint and getGlobalIndexArrayOfBezierPointsInTetrahedron give the same answer
				if(globalIndex!=indexArray[elem]) {
					ADD_FAILURE() << "wrong global index given by  getGlobalIndexOfControlPoint(). It is : "<<globalIndex <<" and should be "<<indexArray[elem]  << std::endl;
					return false;
				}
 				TetrahedronBezierIndex tbi=container->getTetrahedronIndex(elem);
				if(elem!=container->getLocalIndexFromTetrahedronIndex(tbi)) {
					ADD_FAILURE() << "wrong local index given by  getLocalIndexFromTetrahedronIndex(). It is : "<<container->getLocalIndexFromTetrahedronIndex(tbi) <<" and should be "<<elem  << std::endl;
					return false;
				}
				// check that getTetrahedronBezierIndex is consistant with getTetrahedronBezierIndexArray
				if ((tbiArray[elem][0]!=tbi[0]) || (tbiArray[elem][1]!=tbi[1]) || (tbiArray[elem][2]!=tbi[2]) || (tbiArray[elem][3]!=tbi[3])) {
					ADD_FAILURE() << "non consistent indices between getTetrahedronBezierIndexArray() and getTetrahedronBezierIndex(). Got  : "<<tbiArray[elem] <<" versus  "<<tbi  << std::endl;
					return false;
				}
				// check that getLocationFromGlobalIndex is consistent with 
				container->getLocationFromGlobalIndex(globalIndex,location,elementIndex,elementOffset);
				if (elem<4) {
					if ((location!=HighOrderTetrahedronSetTopologyContainer::POINT) || (elementIndex!=container->getTetrahedron(nTetras)[elem]) || (elementOffset!=0)) {
						ADD_FAILURE() << "non consistent indices given by  getLocationFromGlobalIndex() for global index : "<<globalIndex <<std::endl;
						return false;
					}
				}
				else if (elem<(size_t)(4+6*(degree-1))){
					if ((location!=HighOrderTetrahedronSetTopologyContainer::EDGE) || (elementIndex!=container->getEdgesInTetrahedron(nTetras)[(elem-4)/(degree-1)])) {
						ADD_FAILURE() << "non consistent indices given by  getLocationFromGlobalIndex() for global index : "<<globalIndex <<std::endl;
						return false;
					}

				}
				else if (elem<(size_t)(4+6*(degree-1)+2*(degree-1)*(degree-2))){
					size_t nbPointPerEdge=(degree-1)*(degree-2)/2;
					size_t val=(elem-4-6*(degree-1))/(nbPointPerEdge);
					if ((location!=HighOrderTetrahedronSetTopologyContainer::TRIANGLE) || (elementIndex!=container->getTrianglesInTetrahedron(nTetras)[val])) {
						ADD_FAILURE() << "non consistent indices given by  getLocationFromGlobalIndex() for global index : "<<globalIndex <<std::endl;
						return false;
					}
				}
			}

		}
		return( true);
	}
	bool testBezierTetrahedronGeometry()
	{
		HighOrderTetrahedronSetTopologyContainer *container=root->get<HighOrderTetrahedronSetTopologyContainer>(root->SearchDown);
		typename MechanicalObject::SPtr dofs = root->get<MechanicalObject>(std::string("BezierTetrahedronTopology/"));
		typename MechanicalObject::WriteVecCoord coords = dofs->writePositions();
		size_t i,j;	
		BezierDegreeType degree=container->getDegree();

		sofa::helper::vector<TetrahedronBezierIndex> tbiArray=container->getTetrahedronIndexArray();
		for ( i = 0; i<container->getNumberOfTetrahedra(); i++)
		{
			
			const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=container->getGlobalIndexArrayOfControlPoints(i);
			

			for (j=0;j<tbiArray.size();++j) {

				if (j>=4) {
					// test if the position is correct
					Coord pos=coords[indexArray[0]]*(Real)tbiArray[j][0]/degree+coords[indexArray[1]]*(Real)tbiArray[j][1]/degree+coords[indexArray[2]]*(Real)tbiArray[j][2]/degree+coords[indexArray[3]]*(Real)tbiArray[j][3]/degree;
					if ((pos-coords[indexArray[j]]).norm()>1e-5) {
						ADD_FAILURE() << "Wrong control point position in tetrahedron no  : "<<i <<" for point of local index " <<j
						<< " Got point position="<<coords[indexArray[j]]<<" instead of "<<pos<<std::endl;
						return false;
					}
				}

			}
		}
		return true;
	}
	/*
	bool testBezierTetrahedronMass()
	{
		HighOrderTetrahedronSetTopologyContainer *container=root->get<HighOrderTetrahedronSetTopologyContainer>(root->SearchDown);
		BezierTetrahedronSetGeometryAlgorithms *geo=root->get<BezierTetrahedronSetGeometryAlgorithms>(root->SearchDown);
	
		typename MechanicalObject::SPtr dofs = root->get<MechanicalObject>(std::string("BezierTetrahedronTopology/"));
		typename MechanicalObject::WriteVecCoord coords = dofs->writePositions();
		MeshMatrixMass *mass=root->get<MeshMatrixMass>(root->SearchDown);
		const sofa::helper::vector<typename MeshMatrixMass::MassVector> & mv=mass->tetrahedronMassInfo.getValue();
		const sofa::helper::vector<typename MeshMatrixMass::MassType> &ma =mass->vertexMassInfo.getValue();

		size_t i,j,k,rank;	
		BezierDegreeType degree=container->getDegree();
		Real tetraVol1,tetraVol2,totalVol1,totalVol2;
		
		sofa::helper::vector<TetrahedronBezierIndex> tbiArray=container->getTetrahedronBezierIndexArray();
		size_t nbControlPoints=(degree+1)*(degree+2)*(degree+3)/6;
		totalVol1=0;
		for ( i = 0; i<container->getNumberOfTetrahedra(); i++)
		{
			
//				const HighOrderTetrahedronSetTopologyContainer::VecPointID &indexArray=container->getGlobalIndexArrayOfBezierPoints(i);
			
			
			/// get the volume of the tetrahedron
			tetraVol1=geo->computeTetrahedronVolume(i);
			tetraVol2=0;
			// compute the total volume
			totalVol1+=tetraVol1;
			/// check that the sum of the volume matrix elements is equal to the volume of the tetrahedron
			for (rank=0,j=0;j<nbControlPoints;j++) {
				for (k=j;k<nbControlPoints;k++,rank++) {
					if (k==j) 
						// add diagonal term
						tetraVol2+=mv[i][rank]; 
					else 
						// add 2 times off-diagonal term
						tetraVol2+=2*mv[i][rank]; 
				}
			}
			if (fabs(tetraVol1-tetraVol2)>1e-5) {
				ADD_FAILURE() << "Wrong mass matrix in tetrahedron no  : "<<i
				<< " Got total mass="<<tetraVol2<<" instead of "<<tetraVol1<<std::endl;
				return false;
			}
		}
		// compute totalVol2 as the total of the lumped volume
		totalVol2=0;
		for ( i = 0; i<ma.size(); i++)
		{
			totalVol2+=ma[i];
		}
		if (fabs(totalVol1-totalVol2)>1e-5) {
			ADD_FAILURE() << "Wrong total vertex mass value."
			 << " Got total vertex mass="<<totalVol2<<" instead of "<<totalVol1<<std::endl;
			return false;
		}
		return true;
	} */
    void TearDown()
    {
        if (root!=NULL)
            sofa::simulation::getSimulation()->unload(root);
//        cerr<<"tearing down"<<endl;
    }

}; 

// Define the list of DataTypes to instanciate
using testing::Types;
typedef Types<
    Vec3Types
> DataTypes; // the types to instanciate.

// Test suite for all the instanciations
TYPED_TEST_CASE(BezierTetrahedronTopology_test, DataTypes);

// first test topology
TYPED_TEST( BezierTetrahedronTopology_test , testTopology )
{
//	this->loadScene( "tests/SofaTest/BezierTetrahedronTopology.scn");
	this->createScene();
	ASSERT_TRUE( this->testBezierTetrahedronTopology());
	ASSERT_TRUE( this->testBezierTetrahedronGeometry());
//	ASSERT_TRUE( this->testBezierTetrahedronMass());
}



} // namespace sofa
