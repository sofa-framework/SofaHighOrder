#ifndef SOFA_HIGHORDERTOPOLOGY_HOMFFILEEXPORTER_H
#define SOFA_HIGHORDERTOPOLOGY_HOMFFILEEXPORTER_H


#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/objectmodel/DataFileName.h>
#include "HighOrderTetrahedronSetGeometryAlgorithms.h"
#include "HighOrderTriangleSetGeometryAlgorithms.h"
#include <SofaBaseMechanics/MechanicalObject.h>

#include <fstream>
#include "initHighOrderTopology.h"

namespace sofa
{

namespace component
{

namespace misc
{
 using namespace sofa::component::topology;

 template < class DataTypes >
class SOFA_HIGHORDER_TOPOLOGY_API HOMFFileExporter : public core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(HOMFFileExporter, DataTypes),core::objectmodel::BaseObject);

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef sofa::helper::vector<SReal> SeqWeights;
    // the type of shape functions considered
    enum ShapeFunctionType { BEZIER_SHAPE_FUNCTION = 0, LAGRANGE_SHAPE_FUNCTION = 1 };
private:
    sofa::core::topology::BaseMeshTopology* topology;
    sofa::component::container::MechanicalObject<DataTypes>* mstate;
    size_t stepCounter;

    std::ofstream* outfile;

    void writeHOMF();

    HighOrderTetrahedronSetGeometryAlgorithms<DataTypes> *hotsga;
    HighOrderTriangleSetGeometryAlgorithms<DataTypes> *hotrsga;
    HighOrderTriangleSetTopologyContainer *hotrstc;
    HighOrderTetrahedronSetTopologyContainer *hotstc;

public:
    sofa::core::objectmodel::DataFileName filename;
    Data<size_t> formatVersion;	//1 for Simple Legacy Formats
    Data<VecCoord> position;

    Data<unsigned int> exportEveryNbSteps;
    Data<bool> exportAtBegin;
    Data<bool> exportAtEnd;
    Data<bool> overwrite;

    int nbFiles;
    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const HOMFFileExporter<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }
protected:
    HOMFFileExporter();
    virtual ~HOMFFileExporter();
public:
    void init();
    void cleanup();
    void bwdInit();

    void handleEvent(sofa::core::objectmodel::Event *);
};

}

}

}

#endif /* SOFA_HIGHORDERTOPOLOGY_HOMFFILEEXPORTER_H */
