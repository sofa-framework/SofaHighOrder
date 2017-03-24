
#ifndef SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERTETRAHEDRONCONTROLLER_INL
#define SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERTETRAHEDRONCONTROLLER_INL

#include "AdaptiveBezierTetrahedronController.h" 
#include "AdaptiveBezierTetrahedronSetTopologyAlgorithms.h" 
#include "AdaptiveBezierTetrahedronSetTopologyContainer.h" 
#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/AnimateEndEvent.h>
#include <sofa/simulation/Simulation.h>
#include <sofa/simulation/MechanicalVisitor.h>
#include <sofa/simulation/VisualVisitor.h>
#include <sofa/simulation/UpdateMappingVisitor.h>
#include <sofa/defaulttype/DataTypeInfo.h>
#include <sofa/core/objectmodel/Context.h>
#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <sofa/core/objectmodel/KeyreleasedEvent.h>
#include <sofa/core/objectmodel/Data.h>
#include <fstream>
#include <sofa/defaulttype/Vec.h>
#include <cmath>




namespace sofa
{

namespace component
{

namespace controller
{
	typedef 	sofa::component::topology::AdaptiveHighOrderTetrahedronSetTopologyContainer::SeqBezierDegree SeqBezierDegree;

///////////////////////////// Monitor /////////////////////////////////////
template <class DataTypes>
AdaptiveBezierTetrahedronController<DataTypes>::AdaptiveBezierTetrahedronController()
    : d_updateFrequency ( initData ( &d_updateFrequency, "updateFrequency", "number of time step between which the degree of Bezier tetrahedra are modified" ) )
    , d_updateMethod ( initData ( &d_updateMethod,std::string("random"), "updateMethod", "method used to control the degree of Bezier tetrahedra : either \"random\" , \"strain\" or \"stress\" " ) )
    , d_thresholdRaising ( initData ( &d_thresholdRaising, "raisingDegreeValue", "threshold used to control the degree raising of Bezier tetrahedra " ) )
    , d_thresholdLowering ( initData ( &d_thresholdLowering, "loweringDegreeValue", "threshold used to control the degree lowering of Bezier tetrahedra " ) )
{
    if (!f_listening.isSet()) f_listening.setValue(true);


}
/////////////////////////// end Monitor ///////////////////////////////////



////////////////////////////// ~Monitor ///////////////////////////////////
template <class DataTypes>
AdaptiveBezierTetrahedronController<DataTypes>::~AdaptiveBezierTetrahedronController()
{

}
///////////////////////////// end~Monitor /////////////////////////////////




////////////////////////////// init () ////////////////////////////////////
template<class DataTypes>
void AdaptiveBezierTetrahedronController<DataTypes>::init()
{
   if (d_updateMethod==std::string("random")) 
	   updateMethod=RANDOM;

    iterationNumber=0;
//	 std::srand(20);
	getContext()->get(algorithm);
	if (!algorithm) {
		serr << "could not find a AdaptiveBezierTetrahedronSetTopologyAlgorithms object"<<sendl;
	}
	getContext()->get(adaptiveContainer);
	if (!adaptiveContainer) {
		serr << "could not find a AdaptiveBezierTetrahedronSetTopologyContainer object"<<sendl;
	}
	this->f_listening.setValue(true);
	reinit();
}
///////////////////////////// end init () /////////////////////////////////



///////////////////////////// reset () ////////////////////////////////////
template<class DataTypes>
void AdaptiveBezierTetrahedronController<DataTypes>::reset()
{
	 iterationNumber=0;
}
//////////////////////////// end reset () /////////////////////////////////


typedef sofa::core::topology::Topology::TetraID TetraID;

//////////////////////////// reinit () ////////////////////////////////////
template<class DataTypes>
void AdaptiveBezierTetrahedronController<DataTypes>::reinit()
{
   iterationNumber=0;
   // all tetrahedra should be set as linear tetrahedra to start with
  

}
/////////////////////////// end reinit () /////////////////////////////////






template<class DataTypes>
void AdaptiveBezierTetrahedronController<DataTypes>::randomChangeOfDegree() 
{

}
template<class DataTypes>
void AdaptiveBezierTetrahedronController<DataTypes>:: onEndAnimationStep(const double dt) 
{

	if (iterationNumber%(d_updateFrequency.getValue())==0) {


		if (updateMethod==RANDOM) {
			randomChangeOfDegree();
		}
	}
	iterationNumber++;
}
template<class DataTypes>
void AdaptiveBezierTetrahedronController<DataTypes>::onKeyPressedEvent( core::objectmodel::KeypressedEvent* kev )
{
	sofa::helper::vector<TetraID> loweringDegreeTetrahedra,raisingDegreeTetrahedra;
	const  SeqBezierDegree &edgeDegreeArray=adaptiveContainer->getTetrahedronDegreeArray();
	switch(kev->getKey())
	{
		// change the degree of all tetrahedrat
	case 'x':
	case 'X':

		for (size_t i=0;i<adaptiveContainer->getNbTetrahedra();++i) {
			if (edgeDegreeArray[i]==1) 
				raisingDegreeTetrahedra.push_back(i);
			else 
				loweringDegreeTetrahedra.push_back(i);
		}


		algorithm->updateTetrahedronDegree(loweringDegreeTetrahedra,raisingDegreeTetrahedra);

		sofa::simulation::MechanicalPropagatePositionAndVelocityVisitor(sofa::core::MechanicalParams::defaultInstance()).execute(this->getContext());
				sofa::simulation::VisualUpdateVisitor(sofa::core::ExecParams::defaultInstance()).execute(this->getContext());
		sofa::simulation::UpdateMappingVisitor(sofa::core::ExecParams::defaultInstance()).execute(this->getContext());
		break;
	case 'd':
	case 'D':
		{
			// randomly change status of 10% of tetrahedra

			
			size_t nbTetras=adaptiveContainer->getNbTetrahedra();
			size_t nbChanges=nbTetras/10;
			size_t ind;
			std::set<size_t> tetraChangedSet;
			do {
				int random_variable = std::rand();
				ind=nbTetras*random_variable/RAND_MAX;
				// test if index already in set
				if (tetraChangedSet.find(ind)==tetraChangedSet.end()) {
					// if not change its degree
					tetraChangedSet.insert(ind);
					if (edgeDegreeArray[ind]==1) 
						raisingDegreeTetrahedra.push_back(ind);
					else 
						loweringDegreeTetrahedra.push_back(ind);
				}
			} while (tetraChangedSet.size()<nbChanges);

			algorithm->updateTetrahedronDegree(loweringDegreeTetrahedra,raisingDegreeTetrahedra);

			sofa::simulation::MechanicalPropagatePositionAndVelocityVisitor(sofa::core::MechanicalParams::defaultInstance()).execute(this->getContext());
			sofa::simulation::VisualUpdateVisitor(sofa::core::ExecParams::defaultInstance()).execute(this->getContext());
			sofa::simulation::UpdateMappingVisitor(sofa::core::ExecParams::defaultInstance()).execute(this->getContext());
			break;
		}
		
	}

}


/////////////////////////// draw () ////////////////////////////////////
template<class DataTypes>
void AdaptiveBezierTetrahedronController<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
   
}
/////////////////////////// end draw () ////////////////////////////////






} // namespace misc

} // namespace component

} // namespace sofa

#endif
