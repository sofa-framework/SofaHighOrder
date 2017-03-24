#ifndef SOFA_HIGHORDERTOPOLOGY_GENERATEBEZIERCYLINDER_H
#define SOFA_HIGHORDERTOPOLOGY_GENERATEBEZIERCYLINDER_H
#include "initHighOrderTopology.h"

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/topology/BaseMeshTopology.h>

#include <sofa/defaulttype/Vec3Types.h>

namespace sofa
{

namespace component
{

namespace engine
{

/**
 * This class dilates the positions of one DataFields into new positions after applying a dilateation
This dilateation can be either translation, rotation, scale
 */
template <class DataTypes>
class GenerateBezierCylinder : public core::DataEngine
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(GenerateBezierCylinder,DataTypes),core::DataEngine);
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Real Real;
    typedef sofa::core::topology::BaseMeshTopology::SeqTetrahedra SeqTetrahedra;
    typedef typename sofa::core::topology::Topology::Edge Edge;
    typedef sofa::core::topology::BaseMeshTopology::SeqTriangles SeqTriangles;
    typedef typename SeqTetrahedra::value_type Tetrahedron;
    typedef typename SeqTriangles::value_type Triangle;

public:

    GenerateBezierCylinder();

    ~GenerateBezierCylinder() {}

    void init();

    void reinit();

    void update();

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const GenerateBezierCylinder<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

public:
    Data<VecCoord> f_outputTetrahedraPositions; ///< ouput tetrahedra position
	Data<VecCoord> f_outputTrianglesPositions; ///< ouput triangle positions
    Data<SeqTetrahedra> f_tetrahedra; ///< output tetrahedra
	Data<SeqTriangles> f_triangles; ///< output triangles
	Data<sofa::helper::vector<Real> > f_bezierTriangleWeight; ///  output weight for rational Bezier triangles
	Data<sofa::helper::vector<bool> > f_isBezierTriangleRational; ///  for each Bezier triangle indicates if it is rational or integral
    Data<size_t> f_bezierTriangleDegree; /// degree of Bezier triangles
	Data<sofa::helper::vector<Real> > f_bezierTetrahedronWeight; ///  output weight for rational Bezier triangles
    Data<sofa::helper::vector<bool> > f_isBezierTetrahedronRational; ///  for each Bezier tetrahedron indicates if it is rational
	Data<size_t> f_bezierTetrahedronDegree; /// degree of Bezier tetrahedron
    Data<Real > f_radius; /// radius of cylinder 
	Data<Real > f_height; /// height of cylinder
    Data<Coord> f_origin; /// origin
    Data<bool> f_openSurface; /// if the triangulated surface is open or not
    Data<size_t> f_resolutionCircumferential; /// number of points in the circumferential direction
    Data<size_t> f_resolutionRadial; /// number of points in the radial  direction
    Data<size_t> f_resolutionHeight; /// number of points in the height direction
};


#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_ENGINE_GENERATECYLINDER_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_HIGHORDER_TOPOLOGY_API GenerateBezierCylinder<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_HIGHORDER_TOPOLOGY_API GenerateBezierCylinder<defaulttype::Vec3fTypes>;
#endif
#endif

} // namespace engine

} // namespace component

} // namespace sofa

#endif
