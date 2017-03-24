
#ifndef SOFA_COMPONENT_HIGHORDER_TOPOLOGY_INIT_H
#define SOFA_COMPONENT_HIGHORDER_TOPOLOGY_INIT_H

#include <sofa/helper/system/config.h>

#ifdef SOFA_BUILD_SOFAHIGHORDERTOPOLOGY
#  define SOFA_HIGHORDER_TOPOLOGY_API SOFA_EXPORT_DYNAMIC_LIBRARY
#else
#  define SOFA_HIGHORDER_TOPOLOGY_API SOFA_IMPORT_DYNAMIC_LIBRARY
#endif


namespace sofa
{

namespace component
{


void SOFA_HIGHORDER_TOPOLOGY_API initHighOrderTopology();

} // namespace component

} // namespace sofa

#endif

