#ifndef SOFA_COMPONENT_FORCEFIELD_ADAPTIVEBEZIERTETRAHEDRALCOROTATIONALFEMFORCEFIELD_H
#define SOFA_COMPONENT_FORCEFIELD_ADAPTIVEBEZIERTETRAHEDRALCOROTATIONALFEMFORCEFIELD_H

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

#include "HighOrderTetrahedralCorotationalFEMForceField.h"
#include <sofa/core/DataEngine.h>

namespace sofa
{

namespace component
{
namespace topology 
{
	class AdaptiveHighOrderTetrahedronSetTopologyContainer;
}
namespace forcefield
{

using namespace sofa::defaulttype;
using namespace sofa::component::topology;



template<class DataTypes>
class  AdaptiveBezierTetrahedralCorotationalFEMForceField : public HighOrderTetrahedralCorotationalFEMForceField<DataTypes> ,  public core::DataEngine
{
public:
    SOFA_CLASS2(SOFA_TEMPLATE(AdaptiveBezierTetrahedralCorotationalFEMForceField,DataTypes), SOFA_TEMPLATE(HighOrderTetrahedralCorotationalFEMForceField,DataTypes),core::DataEngine);

   
    typedef typename DataTypes::Real        Real        ;
    typedef typename DataTypes::Coord       Coord       ;
    typedef typename DataTypes::Deriv       Deriv       ;
    typedef typename DataTypes::VecCoord    VecCoord    ;
    typedef typename DataTypes::VecDeriv    VecDeriv    ;
    typedef typename DataTypes::VecReal     VecReal     ;
    typedef Data<VecCoord>                  DataVecCoord;
    typedef Data<VecDeriv>                  DataVecDeriv;    



    typedef Mat<3,3,Real>       Mat3x3  ;
    // In case of non 3D template
    typedef Vec<3,Real> Vec3;
    typedef StdVectorTypes< Vec3, Vec3, Real >     GeometricalTypes ; /// assumes the geometry object type is 3D
	typedef typename sofa::component::topology::HighOrderTetrahedronSetGeometryAlgorithms<GeometricalTypes>::VecPointID VecPointID;

	typedef enum { MAX_STRESS, MAX_STRAIN, VON_MISES_STRESS } DeformationIndex;

protected: 
	sofa::component::topology::AdaptiveHighOrderTetrahedronSetTopologyContainer *adaptiveContainer;

   
    AdaptiveBezierTetrahedralCorotationalFEMForceField();

    virtual ~AdaptiveBezierTetrahedralCorotationalFEMForceField();

public:

	virtual void init();
	virtual void handleTopologyChange(core::topology::Topology *topo);
//    void draw(const core::visual::VisualParams* vparams);
	 Data<std::string > f_indexType; /// the string indicating the type of deformation index
     Data<sofa::helper::vector< Real > > f_deformationIndexArray; /// a real array indicating the amount of deformation
protected :
   
	virtual const helper::vector<Mat3x3> &getStiffnessArray(const size_t i,
        typename HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::TetrahedronRestInformation *restTetra);

	virtual const helper::vector<Mat3x3> &getRotatedStiffnessArray(const size_t i,
        const typename HighOrderTetrahedralCorotationalFEMForceField<DataTypes>::TetrahedronRestInformation *restTetra);
public :
	 virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

	 static std::string templateName(const AdaptiveBezierTetrahedralCorotationalFEMForceField<DataTypes>* = NULL)
	 {
		 return DataTypes::Name();
	 }
	  template< class T>
    static std::string shortName( const T* /*ptr*/ = NULL, core::objectmodel::BaseObjectDescription* = NULL )
    {
        return std::string("AdaptiveBezierTetrahedralCorotationalFEMForceField");
    }
	 virtual void update(); // callback when the engine is used
};


} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_FASTTETRAHEDRALCOROTATIONALFORCEFIELD_H
