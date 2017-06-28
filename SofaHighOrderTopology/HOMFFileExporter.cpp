#include "HOMFFileExporter.inl"

#include "initHighOrderTopology.h"
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/ObjectFactory.h>
namespace sofa
{

namespace component
{

namespace misc
{
using namespace sofa::defaulttype;

SOFA_DECL_CLASS(HOMFFileExporter)

int HOMFFileExporterClass = core::RegisterObject("Write High Order Mesh in HOMF format")
#ifdef SOFA_FLOAT
    .add< HOMFFileExporter<Vec3fTypes> >(true)
    .add< HOMFFileExporter<Vec2fTypes> >()
#else
    .add< HOMFFileExporter<Vec3dTypes> >(true)
    .add< HOMFFileExporter<Vec2dTypes> >()
#ifndef SOFA_DOUBLE
    .add< HOMFFileExporter<Vec3fTypes> >() // default template
    .add< HOMFFileExporter<Vec2fTypes> >()
#endif
#endif
;

#ifndef SOFA_FLOAT
template class SOFA_HIGHORDER_TOPOLOGY_API HOMFFileExporter<Vec3dTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API HOMFFileExporter<Vec2dTypes>;
#endif

#ifndef SOFA_DOUBLE
template class SOFA_HIGHORDER_TOPOLOGY_API HOMFFileExporter<Vec3fTypes>;
template class SOFA_HIGHORDER_TOPOLOGY_API HOMFFileExporter<Vec2fTypes>;
#endif


}

}

}
