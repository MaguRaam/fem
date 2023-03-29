#include <fstream>
#include <stdexcept>
#include <iomanip>
#include "mesh.h"

// read coordinates, connectivity and boundary nodes from file:
void ReadMesh(std::string coord_filename, std::string conn_filename, std::string bnode_filename, Mesh &M)
{
    // read coordinates:
    std::ifstream file(coord_filename);
    if (!file.good())
        throw std::runtime_error("Failed to open " + coord_filename);

    // get no of nodes:
    file >> M.nNodes;

    // reserve coordinate vector:
    M.coordinates.reserve(2 * M.nNodes);

    // loop over nodes and push back x, y:
    double x, y;

    for (int n = 0; n < M.nNodes; ++n)
    {
        file >> x >> y;
        M.coordinates.push_back(x);
        M.coordinates.push_back(y);
    }

    file.close();

    // read connectivity:
    file.open(conn_filename);
    if (!file.good())
        throw std::runtime_error("Failed to open " + conn_filename);

    // get no of elements:
    file >> M.nElements;

    // reserve connectivity vector:
    M.connectivity.reserve(3 * M.nElements);

    // loop over elements and push back node1, node2, node3;
    int node1, node2, node3;

    for (int e = 0; e < M.nElements; ++e)
    {
        file >> node1 >> node2 >> node3;
        M.connectivity.push_back(node1);
        M.connectivity.push_back(node2);
        M.connectivity.push_back(node3);
    }

    file.close();

    // read boundary nodes
    file.open(bnode_filename);
    if (!file.good())
        throw std::runtime_error("Failed to open " + bnode_filename);

    // get no of boundary nodes:
    file >> M.nBnodes;

    // reserve boundary nodes vector:
    M.bdNodes.reserve(M.nBnodes);

    int node;
    for (int n = 0; n < M.nBnodes; ++n)
    {
        file >> node;
        M.bdNodes.push_back(node);
    }

    file.close();
}

// give local to global map:  returns the global index corresponding to a given node of an element
int Local2GlobalMap(const Mesh &M, int elm_num, int loc_node_num)
{
    return M.connectivity[3 * elm_num + loc_node_num];
}

// write solution in a format that can be visualized using python:
void WriteSolution(std::string filename, const Mesh &M, const double *U)
{
    std::ofstream file(filename);
    if (!file.good())
        throw std::runtime_error("Failed to open " + filename);

    // write header:
    file << M.nNodes << " " << M.nElements << "\n";

    // write nodal coordinates:
    for (int n = 0; n < M.nNodes; ++n)
        file << M.coordinates[2 * n] << " " << M.coordinates[2 * n + 1] << "\n";

    // write connectivity:
    for (int e = 0; e < M.nElements; ++e)
        file << M.connectivity[3 * e] << " " << M.connectivity[3 * e + 1] << " " << M.connectivity[3 * e + 2] << "\n";

    // write solution at nodes:
    for (int n = 0; n < M.nNodes; ++n)
        file << std::scientific << std::setprecision(6) << U[n] << "\n";
}