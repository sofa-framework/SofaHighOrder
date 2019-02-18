#include <sofa/core/ObjectFactory.h>
#include "MeshHOMFLoader.h"
#include <sofa/core/visual/VisualParams.h>
#include <iostream>
#include <fstream>

namespace sofa
{

namespace component
{

namespace loader
{

using namespace sofa::defaulttype;
using std::string;
using std::stringstream;

SOFA_DECL_CLASS(MeshHOMFLoader)

int MeshHOMFLoaderClass = core::RegisterObject("Specific mesh loader for HOMF file format.")
        .add< MeshHOMFLoader >()
        ;

MeshHOMFLoader::MeshHOMFLoader() : MeshLoader()
, d_positions2D(initData(&d_positions2D, "position2D", "Array of 2d positions of control points"))
, d_weights(initData(&d_weights, "weights", "Array of weights for rational splines"))
, d_degree(initData(&d_degree, "degree", "degree of spline mesh"))
, d_swapEdges(initData(&d_swapEdges, (bool)false,"swapEdges", "swap order of edges"))

{
    d_positions2D.setPersistent(false);
    d_weights.setPersistent(false);
    d_degree.setPersistent(false);

}

bool MeshHOMFLoader::load()
{
    sout << "Loading HOMF file: " << m_filename << sendl;

    string cmd;
    bool fileRead = false;
    unsigned int homfFormat = 0;

    // -- Loading file
    const char* filename = m_filename.getFullPath().c_str();
    std::ifstream file(filename);

    if (!canLoad())
        return false;

    // -- Looking for HOMF version of this file.
    std::getline(file, cmd); //Version

    if (cmd.compare(0, 12, "HOMF Version") == 0) {
        string version = cmd.substr(13, 1);
        homfFormat = std::stoi(version);
        if (homfFormat > 1) {
            serr << "Wrong version number "<< homfFormat<< ". Closing File" << sendl;
            file.close();
            return false;
        }
        fileRead = readHOMF(file, homfFormat);
        return fileRead;
    }
    else {
        serr << "Wrong format. Closing File" << sendl;
        file.close();
        return false;
    }
   
}

bool MeshHOMFLoader::readHOMF(std::ifstream &file, const unsigned int homfFormat)
{
    sout << "Reading HOMF file: " << homfFormat << sendl;

    string cmd;

    // ---- Loading header - part 1 ----
    // read dimension of embedding space (either 2 for 2d positions or 3 for 3d positions)
    size_t dimensionEmbedding = 0;
    // read dimension of simplex (either 2 for triangles or 3 for tetrahedra)
    size_t dimensionSimplex = 0;
    file >> dimensionEmbedding >> dimensionSimplex;
    if ((dimensionEmbedding != 2) && (dimensionEmbedding != 3)) {
        serr << "Wrong embedding dimension  " << dimensionEmbedding << ". Closing File" << sendl;
        file.close();
        return false;
    }
    if ((dimensionSimplex != 2) && (dimensionSimplex != 3)) {
        serr << "Wrong simplex dimension  " << dimensionSimplex << ". Closing File" << sendl;
        file.close();
        return false;
    }
    // ---- Loading header - part 2 ----
    size_t my_degree;
    file >> my_degree;
    d_degree.setValue(my_degree);

    ShapeFunctionType sfType = BEZIER_SHAPE_FUNCTION;
    size_t bsf;
    file >> bsf;
    sfType = (ShapeFunctionType)bsf;

    if (bsf > 1) {
        serr << "Wrong shape function type  " << sfType << ". Closing File" << sendl;
        file.close();
        return false;
    }
    bool rationalMesh = true;
  
    // the number of coordinates to be read
    size_t numberInputCoordinates = dimensionEmbedding + (rationalMesh == true);

    // --- Loading Vertices of macro simplicial mesh ---
    size_t numberVertices = 0;
    file >> numberVertices;

    // read each input coordinate for each vertex;
    double x, y, z, w;

    helper::vector<sofa::defaulttype::Vector3>& my_positions = *(d_positions.beginEdit());
    helper::vector<sofa::defaulttype::Vector2>& my_positions2D = *(d_positions2D.beginEdit());
    helper::vector<SReal>& my_weights = *(d_weights.beginEdit());

    size_t i;
    for (i = 0; i < numberVertices; ++i) {

        switch (numberInputCoordinates) {
        case 2:
            file >> x >> y;
            break;
        case 3:
            file >> x >> y >> z;
            break;
        case 4:
            file >> x >> y >> z >> w;
            break;
        }
        if (dimensionEmbedding == 2) {
            my_positions2D.push_back(Vector2(x, y));
            my_positions.push_back(Vector3(x, y,0));
        }
        else {
            my_positions.push_back(Vector3(x, y, z));
        }
        if (rationalMesh) {
            if (dimensionEmbedding == 2) {
                my_weights.push_back((SReal)z);
            }
            else
                my_weights.push_back((SReal)w);
        }
    }
    d_positions.endEdit();
    d_positions2D.endEdit();
    d_weights.endEdit();


    // --- Loading Simplices of macro simplicial mesh ---

    // ----- load edges ------ 
    helper::vector<Edge>& my_edges = *(d_edges.beginEdit());
    size_t numberEdges = 0;
    file >> numberEdges;

    size_t v0, v1, v2, v3;

    for (i = 0; i < numberEdges; ++i) {
        file >> v0 >> v1;
        addEdge(&my_edges, Edge(v0, v1));
    }
    d_edges.endEdit();

    // ----- load triangles ------ 
    size_t numberTriangles = 0;
    file >> numberTriangles;

    helper::vector<Triangle>& my_triangles = *(d_triangles.beginEdit());

    for (i = 0; i < numberTriangles; ++i) {
        file >> v0 >> v1 >> v2;
        addTriangle(&my_triangles, Triangle(v0, v1, v2));
    }
    d_triangles.endEdit();

    // ----- load tetrahedra ------ 
    size_t numberTetrahedra = 0;
    if (dimensionSimplex == 3) {
       
        file >> numberTetrahedra;

        helper::vector<Tetrahedron>& my_tetrahedra = *(d_tetrahedra.beginEdit());

        for (i = 0; i < numberTetrahedra; ++i) {
            file >> v0 >> v1 >> v2 >> v3;
            addTetrahedron(&my_tetrahedra, Tetrahedron(v0, v1, v2, v3));
        }
        d_tetrahedra.endEdit();
    }


    // --- Loading index position of control points  ---
    size_t i0, i1, i2, i3, i4, i5;
    /// loading edge control points indices
    if (my_degree > 1) {
        helper::vector<HighOrderEdgePosition >& my_highOrderEdgePositions = *(d_highOrderEdgePositions.beginEdit());
        size_t numberEdgeControlPoints = (my_degree - 1)*numberEdges;


        for (i = 0; i < numberEdgeControlPoints; ++i) {
            file >> i0 >> i1 >> i2 >> i3;
			if (d_swapEdges.getValue())
				my_highOrderEdgePositions.push_back(HighOrderEdgePosition(i0, i1, i3, i2));
			else
	           my_highOrderEdgePositions.push_back(HighOrderEdgePosition(i0, i1, i2, i3));
	
        }
        d_highOrderEdgePositions.endEdit();
    }
    /// loading triangle control points indices
    if (my_degree > 2) {
        helper::vector<HighOrderTrianglePosition >& my_highOrderTrianglePositions = *(d_highOrderTrianglePositions.beginEdit());
        size_t numberTriangleControlPoints = (my_degree - 1)*(my_degree - 2)*numberTriangles/2;


        for (i = 0; i < numberTriangleControlPoints; ++i) {
            file >> i0 >> i1 >> i2 >> i3 >> i4;
            my_highOrderTrianglePositions.push_back(HighOrderTrianglePosition(i0, i1, i2, i3, i4));
        }
        d_highOrderTrianglePositions.endEdit();
    }
    /// loading triangle control points indices
    if ((my_degree > 3) && (dimensionSimplex>2)) {
        helper::vector<HighOrderTetrahedronPosition >& my_highOrderTetrahedronPositions = *(d_highOrderTetrahedronPositions.beginEdit());
        size_t numberTetrahedronControlPoints = (my_degree - 1)*(my_degree - 2)*(my_degree - 3)*numberTetrahedra / 6;


        for (i = 0; i < numberTetrahedronControlPoints; ++i) {
            file >> i0 >> i1 >> i2 >> i3 >> i4 >> i5;
            my_highOrderTetrahedronPositions.push_back(HighOrderTetrahedronPosition(i0, i1, i2, i3, i4, i5));
        }
        d_highOrderTetrahedronPositions.endEdit();
    }
  

    file.close();
    return true;
}


} // namespace loader

} // namespace component

} // namespace sofa

