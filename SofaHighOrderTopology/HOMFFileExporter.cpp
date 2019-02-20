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
    .add< HOMFFileExporter<Vec3Types> >(true)
    .add< HOMFFileExporter<Vec2Types> >()
;

template class SOFA_HIGHORDER_TOPOLOGY_API HOMFFileExporter<Vec3Types>;
template class SOFA_HIGHORDER_TOPOLOGY_API HOMFFileExporter<Vec2Types>;

}

}

}
