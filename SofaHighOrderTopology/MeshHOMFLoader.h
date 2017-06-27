#ifndef SOFA_HIGHORDERTOPOLOGY_MESHHOMFLOADER_H
#define SOFA_HIGHORDERTOPOLOGY_MESHHOMFLOADER_H


#include <sofa/core/loader/MeshLoader.h>
#include "initHighOrderTopology.h"
namespace sofa
{

namespace component
{

namespace loader
{

class SOFA_HIGHORDER_TOPOLOGY_API MeshHOMFLoader : public sofa::core::loader::MeshLoader
{
public:
    SOFA_CLASS(MeshHOMFLoader,sofa::core::loader::MeshLoader);
    // the type of shape functions considered
    enum ShapeFunctionType {BEZIER_SHAPE_FUNCTION=0, LAGRANGE_SHAPE_FUNCTION=1};
    virtual bool load();

    template <class T>
    static bool canCreate ( T*& obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg )
    {
        return BaseLoader::canCreate (obj, context, arg);
    }

    // Point coordinates in 2D in double.
    Data< helper::vector<sofa::defaulttype::Vec<2, SReal> > > d_positions2D;
    // Point weights for rational splines
    Data< helper::vector<SReal>  > d_weights;

    // the degree of the spline
    Data< size_t  > d_degree;

protected:

    bool readHOMF(std::ifstream &file, const unsigned int homfFormat);

    void addPosition2D(helper::vector< sofa::defaulttype::Vec<2, SReal> >* pPositions, const sofa::defaulttype::Vec<2, SReal> &p);
    void addPosition2D(helper::vector<sofa::defaulttype::Vec<2, SReal> >* pPositions, SReal x, SReal y);

    void addWeight(helper::vector<SReal >* pWeights, SReal w);

public:

};




} // namespace loader

} // namespace component

} // namespace sofa

#endif
