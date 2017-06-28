#include "HOMFFileExporter.h"

#include <sstream>



#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/simulation/AnimateEndEvent.h>
#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <sofa/core/objectmodel/KeyreleasedEvent.h>
#include "BezierTetrahedronSetGeometryAlgorithms.h"
#include "BezierTriangleSetGeometryAlgorithms.h"
namespace sofa
{

namespace component
{

namespace misc
{


template< class DataTypes>
HOMFFileExporter< DataTypes >::HOMFFileExporter()
    : stepCounter(0), outfile(NULL)
    , filename( initData(&filename, "filename", "output HOMF file name"))
    , formatVersion( initData(&formatVersion, (size_t) 1, "formatVersion", "Set the the format version"))
    , position( initData(&position, "position", "points position (will use points from topology or mechanical state if this is empty)"))
    , exportEveryNbSteps( initData(&exportEveryNbSteps, (unsigned int)0, "exportEveryNumberOfSteps", "export file only at specified number of steps (0=disable)"))
    , exportAtBegin( initData(&exportAtBegin, false, "exportAtBegin", "export file at the initialization"))
    , exportAtEnd( initData(&exportAtEnd, false, "exportAtEnd", "export file when the simulation is finished"))
    , overwrite( initData(&overwrite, false, "overwrite", "overwrite the file, otherwise create a new file at each export, with suffix in the filename"))
{
}
template< class DataTypes>
HOMFFileExporter< DataTypes >::~HOMFFileExporter()
{
    if (outfile)
        delete outfile;
}
template< class DataTypes>
void HOMFFileExporter< DataTypes >::init()
{
    sofa::core::objectmodel::BaseContext* context = this->getContext();
    context->get(topology);
    context->get(mstate);

    if (!topology)
    {
        serr << "HOMFFileExporter : error, no topology ." << sendl;
        return;
    }
    else sout << "HOMFFileExporter: found topology " << topology->getName() << sendl;


 

    hotsga = NULL;
    this->getContext()->get(hotsga, sofa::core::objectmodel::BaseContext::SearchUp);
    hotrsga = NULL;
    this->getContext()->get(hotrsga, sofa::core::objectmodel::BaseContext::SearchUp);

    if ((!hotrsga) && (!hotsga)) {
        serr << " Could not find a HighOrderTetrahedronSetGeometryAlgorithms or HighOrderTriangleSetGeometryAlgorithms object" << sendl;
        exit(1);
    }
    if (hotsga) {
        hotstc = NULL;
        this->getContext()->get(hotstc, sofa::core::objectmodel::BaseContext::SearchUp);
        if ((!hotstc)) {
            serr << " Could not find a HighOrderTetrahedronSetTopologyContainer object" << sendl;
            exit(1);
        }
    }
    else {
        hotrstc = NULL;
        this->getContext()->get(hotrstc, sofa::core::objectmodel::BaseContext::SearchUp);
        if ((!hotrstc)) {
            serr << " Could not find a HighOrderTriangleSetTopologyContainer object" << sendl;
            exit(1);
        }
    }
    nbFiles = 0;
}

template< class DataTypes>
void HOMFFileExporter< DataTypes >::writeHOMF()
{
    std::string _filename = filename.getFullPath();

    std::ostringstream oss;
    oss << "_" << nbFiles;

    if (_filename.size() > 3) {
        std::string ext;
        std::string baseName;
        if (_filename.substr(_filename.size()-4)==".hom") {
            ext = ".hom";
            baseName = _filename.substr(0, _filename.size()-4);
        }



        /// no extension given => default "vtu"
        if (ext == "") {
            ext = ".hom";
            baseName = _filename;
        }

        if (overwrite.getValue())
            _filename = baseName + ext;
        else
            _filename = baseName + oss.str() + ext;

    }


    outfile = new std::ofstream(_filename.c_str());
    if( !outfile->is_open() )
    {
        serr << "Error creating file "<<_filename<<sendl;
        delete outfile;
        outfile = NULL;
        return;
    }
    const typename DataTypes::VecCoord& p = (mstate->read(core::ConstVecCoordId::position())->getValue());
 //  helper::ReadAccessor<Data<VecCoord> > pointsPos = position;

    const size_t nbp = (!p.empty()) ? p.size() : topology->getNbPoints();

    //Write header
    *outfile << "HOMF Version 1" << std::endl;

    //write embedding dimension and simplex dimension
    //write also the degree of the spline and the type of shape function
    bool rationalMesh;
    size_t my_degree;
    if (hotsga) {
        *outfile << Coord::spatial_dimensions << " 3" << std::endl;
        *outfile << (size_t)hotstc->getDegree() << std::endl;
        my_degree = (size_t)hotstc->getDegree();
        if (dynamic_cast<BezierTetrahedronSetGeometryAlgorithms<DataTypes> *>(hotsga) != NULL) {
            *outfile << (size_t)BEZIER_SHAPE_FUNCTION << std::endl;
        }
        else {
            *outfile << (size_t)LAGRANGE_SHAPE_FUNCTION << std::endl;
        }

        rationalMesh = true;
    }
    else {
        *outfile << Coord::spatial_dimensions << " 2" << std::endl;
        *outfile << (size_t)hotrstc->getDegree() << std::endl;
        my_degree = (size_t)hotrstc->getDegree();
        if (dynamic_cast<BezierTriangleSetGeometryAlgorithms<DataTypes> *>(hotrsga) != NULL) {
            *outfile << (size_t)BEZIER_SHAPE_FUNCTION << std::endl;
        }
        else {
            *outfile << (size_t)LAGRANGE_SHAPE_FUNCTION << std::endl;
        }
        rationalMesh = true;
    }


    // write the number of control points
    *outfile << nbp <<  std::endl;

    size_t i;
    if (hotsga) {
        const SeqWeights &weights = hotstc->getWeightArray();

        for (i = 0; i < nbp; ++i) {
            *outfile << p[i] << " " << weights[i] << std::endl;
        }
    }
    else {
        const SeqWeights &weights = hotrstc->getWeightArray();

        for (i = 0; i < nbp; ++i) {
            *outfile << p[i] << " " << weights[i] << std::endl;
        }
    }

    
    // export edges
    const sofa::core::topology::BaseMeshTopology::SeqEdges &edges = topology->getEdges();
    size_t nedges = edges.size();
    *outfile << nedges << std::endl;
    for (i = 0; i < nedges; ++i) {
        *outfile << edges[i][0] << " "<< edges[i][1]<< std::endl;
    }
    // export triangles
    const sofa::core::topology::BaseMeshTopology::SeqTriangles &triangles = topology->getTriangles();
    size_t ntriangles = triangles.size();
    *outfile << ntriangles << std::endl;
    for (i = 0; i < ntriangles; ++i) {
        *outfile << triangles[i][0] << " " << triangles[i][1] << " "<< triangles[i][2]<< std::endl;
    }
    // export tetrahedra
    size_t ntetras=0;
    if (hotsga) {
        const sofa::core::topology::BaseMeshTopology::SeqTetrahedra &tetras = topology->getTetrahedra();
        ntetras = tetras.size();
        *outfile << ntetras << std::endl;
        for (i = 0; i < ntetras; ++i) {
            *outfile << tetras[i][0] << " " << tetras[i][1] << " " << tetras[i][2] << " "<< tetras[i][3]<< std::endl;
        }
    }
    // --- Writing index position of control points  ---
    if ((hotrsga)&&(my_degree>1)) {
        size_t nTriangleVertices= hotrstc->getNumberOfTriangularPoints();
        sofa::component::topology::HighOrderTriangleSetTopologyContainer::HighOrderTrianglePointLocation location;
        size_t elementIndex;
        size_t elementOffset;
        typedef sofa::helper::fixed_array<size_t, 3> elementInfo;
        std::set<elementInfo> edgeInfo;
        std::set<elementInfo> triangleInfo;
        for (i = nTriangleVertices; i < nbp; ++i) {
            hotrstc->getLocationFromGlobalIndex(i, location, elementIndex,elementOffset);
            if (location == HighOrderTriangleSetTopologyContainer::EDGE)
                edgeInfo.insert(elementInfo(i, elementIndex, elementOffset));
            else if (location == HighOrderTriangleSetTopologyContainer::TRIANGLE)
                triangleInfo.insert(elementInfo(i, elementIndex, elementOffset));
        }
        assert(edgeInfo.size() == (my_degree - 1)*nedges);
        assert(triangleInfo.size() == (my_degree - 1)* (my_degree - 2)*ntriangles/2);

        std::set<elementInfo>::iterator ite;
        sofa::component::topology::EdgeIndexVector eiv;
        for (ite = edgeInfo.begin(); ite != edgeInfo.end(); ++ite) {
            hotrstc->getEdgeIndexVectorFromEdgeOffset((*ite)[2], eiv);
            *outfile << (*ite)[0] << " " << (*ite)[1] << " " << (size_t)eiv[0] << " " << (size_t)eiv[1] << std::endl;
        }
        if (my_degree > 2) {
            std::set<elementInfo>::iterator itt;
            sofa::component::topology::TriangleIndexVector tiv;
            for (itt = triangleInfo.begin(); itt != triangleInfo.end(); ++itt) {
                hotrstc->getTriangleIndexVectorFromTriangleOffset((*itt)[2], tiv);
                *outfile << (*itt)[0] << " " << (*itt)[1] << " " <<(size_t) tiv[0] << " " << (size_t)tiv[1] << " "<< (size_t)tiv[2] << std::endl;
            }
        }
    }
    else if ((hotsga) && (my_degree > 1)) {
        size_t nTetrahedronVertices = hotstc->getNumberOfTetrahedralPoints();
        sofa::component::topology::HighOrderTetrahedronSetTopologyContainer::HighOrderTetrahedronPointLocation location;
        size_t elementIndex;
        size_t elementOffset;
        typedef sofa::helper::fixed_array<size_t, 3> elementInfo;
        std::set<elementInfo> edgeInfo;
        std::set<elementInfo> triangleInfo;
        std::set<elementInfo> tetraInfo;
        for (i = nTetrahedronVertices; i < nbp; ++i) {
            hotstc->getLocationFromGlobalIndex(i, location, elementIndex, elementOffset);
            if (location == HighOrderTetrahedronSetTopologyContainer::EDGE)
                edgeInfo.insert(elementInfo(i, elementIndex, elementOffset));
            else if (location == HighOrderTetrahedronSetTopologyContainer::TRIANGLE)
                triangleInfo.insert(elementInfo(i, elementIndex, elementOffset));
            else if (location == HighOrderTetrahedronSetTopologyContainer::TETRAHEDRON)
                tetraInfo.insert(elementInfo(i, elementIndex, elementOffset));
        }
        assert(edgeInfo.size() == (my_degree - 1)*nedges);
        assert(triangleInfo.size() == (my_degree - 1)* (my_degree - 2)*ntriangles / 2);
        assert(tetraInfo.size() == (my_degree - 1)* (my_degree - 2)*(my_degree - 3)*ntetras / 6);

        std::set<elementInfo>::iterator ite;
        sofa::component::topology::EdgeIndexVector eiv;
        for (ite = edgeInfo.begin(); ite != edgeInfo.end(); ++ite) {
            hotstc->getEdgeIndexVectorFromEdgeOffset((*ite)[2], eiv);
            *outfile << (*ite)[0] << " " << (*ite)[1] << " " << (size_t)eiv[0] << " " << (size_t)eiv[1] << std::endl;
        }
        if (my_degree > 2) {
            std::set<elementInfo>::iterator itt;
            sofa::component::topology::TriangleIndexVector tiv;
            for (itt = triangleInfo.begin(); itt != triangleInfo.end(); ++itt) {
                hotstc->getTriangleIndexVectorFromTriangleOffset((*itt)[2], tiv);
                *outfile << (*itt)[0] << " " << (*itt)[1] << " " << (size_t)tiv[0] << " " << (size_t)tiv[1] << " " << (size_t)tiv[2] << std::endl;
            }
        }
        if (my_degree > 3) {
            std::set<elementInfo>::iterator ittet;
            sofa::component::topology::TetrahedronIndexVector ttiv;
            for (ittet = tetraInfo.begin(); ittet != tetraInfo.end(); ++ittet) {
                hotstc->getTetrahedronIndexVectorFromTetrahedronOffset((*ittet)[2], ttiv);
                *outfile << (*ittet)[0] << " " << (*ittet)[1] << " " << (size_t)ttiv[0] << " " << (size_t)ttiv[1] << " " << (size_t)ttiv[2] << " " << (size_t)ttiv[3] << std::endl;
            }
        }
    }

    outfile->close();
    sout << filename << " written" << sendl;

    ++nbFiles;
}

template< class DataTypes>
void HOMFFileExporter< DataTypes >::handleEvent(sofa::core::objectmodel::Event *event)
{
    if (sofa::core::objectmodel::KeypressedEvent::checkEventType(event))
    {
        sofa::core::objectmodel::KeypressedEvent* ev = static_cast<sofa::core::objectmodel::KeypressedEvent*>(event);

        std::cout << "key pressed " << std::endl;
        switch(ev->getKey())
        {

        case 'W':
        case 'w':
            writeHOMF();
            break;

        }
    }


    if ( /*simulation::AnimateEndEvent* ev =*/ simulation::AnimateEndEvent::checkEventType(event))
    {
        unsigned int maxStep = exportEveryNbSteps.getValue();
        if (maxStep == 0) return;

        stepCounter++;
        if(stepCounter >= maxStep)
        {
            stepCounter = 0;
            writeHOMF();
        }
    }
}
template< class DataTypes>
void HOMFFileExporter< DataTypes >::cleanup()
{
    if (exportAtEnd.getValue())
        writeHOMF();

}
template< class DataTypes>
void HOMFFileExporter< DataTypes >::bwdInit()
{
    if (exportAtBegin.getValue())
        writeHOMF();
}

}

}

}
