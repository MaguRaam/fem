// 2d L2 projection
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <functional>
#include <tuple>
#include <cassert>
#include <Eigen/Sparse>
#include <Eigen/Dense>

struct point
{
    double x, y;
};

void read_mesh(const std::string &coord_filename, std::vector<point> &coordinates, const std::string &conn_filename, std::vector<std::array<int, 3>> &connectivity)
{

    // open coordinates file:
    std::ifstream infile(coord_filename, std::ios::in);
    assert(infile.is_open() && "ReadMesh:: Coordinates file does not exist");
    coordinates.clear();

    // coordinates:
    double x, y;

    // loop over file line by line and push x and y to coordinates vector:
    while (infile >> x >> y)
        coordinates.push_back({x, y});

    infile.close();

    // open connectivity file:
    infile.open(conn_filename, std::ios::in);
    assert(infile.good() && "ReadMesh:: Connectivity file does not exist");
    connectivity.clear();

    // global node no of triangles:
    int node1, node2, node3;

    while (infile >> node1 >> node2 >> node3)
        connectivity.push_back({node1, node2, node3});

    infile.close();
}

void write_mesh(const std::string &filename, std::vector<point> &coordinates, std::vector<std::array<int, 3>> &connectivity, const double *nodal_field = nullptr)
{

    std::ofstream outfile(filename, std::ios::out);
    outfile.flags(std::ios::dec | std::ios::scientific);
    outfile.precision(16);

    // Line 1:
    if (nodal_field == nullptr)
        outfile << "VARIABLES = \"X\", \"Y\" \n";
    else
        outfile << "VARIABLES = \"X\", \"Y\" \"F1\" \n";

    // Line 2:
    const int nNodes = coordinates.size();
    const int nElements = connectivity.size();
    outfile << "ZONE t=\"t:0\", N=" << nNodes << ", E=" << nElements << ", F=FEPOINT, ET=TRIANGLE\n";

    // Nodal coordinates and field values
    for (int n = 0; n < nNodes; ++n)
    {
        if (nodal_field != nullptr)
            outfile << coordinates[n].x << " " << coordinates[n].y << " " << nodal_field[n] << "\n";
        else
            outfile << coordinates[n].x << " " << coordinates[n].y << "\n";
    }

    // Connectivity
    for (int e = 0; e < nElements; ++e)
        outfile << connectivity[e][0] + 1 << " " << connectivity[e][1] + 1 << " " << connectivity[e][2] + 1 << "\n";

    outfile.close();
}