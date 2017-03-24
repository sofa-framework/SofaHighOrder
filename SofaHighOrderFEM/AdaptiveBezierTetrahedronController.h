#ifndef SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERTETRAHEDRONCONTROLLER_H
#define SOFA_HIGHORDERDFEM_ADAPTIVEBEZIERTETRAHEDRONCONTROLLER_H



#include "initHighOrderFEM.h"
#include <SofaUserInteraction/Controller.h>

namespace sofa
{

namespace component
{
	namespace topology { template<typename  D>  class AdaptiveBezierTetrahedronSetTopologyAlgorithms; } 
	namespace topology {  class AdaptiveHighOrderTetrahedronSetTopologyContainer; } 
	
	

namespace controller
{


template <class DataTypes>
class  SOFA_HIGHORDER_FEM_API  AdaptiveBezierTetrahedronController:  public  Controller
{
	  
public:
	 SOFA_CLASS(SOFA_TEMPLATE(AdaptiveBezierTetrahedronController,DataTypes), Controller);

    typedef typename DataTypes::VecReal VecReal;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;

	typedef enum { RANDOM=1, STRAIN, STRESS} UpdateMethod;
protected:
    AdaptiveBezierTetrahedronController ();
    ~AdaptiveBezierTetrahedronController ();
public:
	/// the frequency at which the order of Bezier tetrahedra are modified
	Data<size_t> d_updateFrequency;
	/// method used to change the degree of Bezier tetrahedra
	Data<std::string> d_updateMethod;
	UpdateMethod updateMethod;
	/// thresholds used to control the degree raising (resp lowering) of Bezier tetrahedron used in various methods for 
	Data<Real> d_thresholdRaising;
	Data<Real> d_thresholdLowering;


public:
    //init data
    virtual void init ();

    //reset Monitored values
    virtual void reset ();

    /**initialize gnuplot files
    *called when ExportGnuplot box is checked
    */
    virtual void reinit();

    /**function called at every step of simulation;
  
    */
    virtual void onEndAnimationStep(const double dt);

	 virtual void onKeyPressedEvent(core::objectmodel::KeypressedEvent *);

    virtual void draw (const core::visual::VisualParams* vparams);

  

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const AdaptiveBezierTetrahedronController<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

	void randomChangeOfDegree();
protected :
	sofa::component::topology::AdaptiveBezierTetrahedronSetTopologyAlgorithms<DataTypes> *algorithm;
	sofa::component::topology::AdaptiveHighOrderTetrahedronSetTopologyContainer *adaptiveContainer;
	size_t iterationNumber;

};


} // namespace misc

} // namespace component

} // namespace sofa

#endif
