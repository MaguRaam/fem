#ifndef MESH_H
#define MESH_H

#include <string>
#include <vector>

// structure to store triangle mesh:
struct Mesh
{
    int nNodes;
    std::vector<double> coordinates;
    int nElements;
    std::vector<int> connectivity;
    int nBnodes;
    std::vector<int> bdNodes;
};

void ReadMesh(std::string, std::string, std::string, Mesh&);
int Local2GlobalMap(const Mesh&, int, int);
void WriteSolution(std::string, const Mesh&, const double*);

#endif
