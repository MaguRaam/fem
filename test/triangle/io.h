#pragma once

#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>

// read mesh:
void read_mesh(const std::string coord_filename,
               std::vector<double> &coordinates,
               const std::string conn_filename,
               std::vector<int> &connectivity)
{
    // open coordinates file:
    std::ifstream infile(coord_filename, std::ios::in);
    assert(infile.is_open() && "ReadMesh:: Coordinates file does not exist");
    coordinates.clear();

    // coordinates:
    double x, y;

    // loop over file line by line and push x and y to coordinates vector:
    while (infile >> x >> y)
    {
        coordinates.push_back(x);
        coordinates.push_back(y);
    }

    infile.close();

    // compute no of nodes:
    assert(coordinates.size() % 2 == 0 && "x,y pair doesn't exist");
    const int nNodes = static_cast<int>(coordinates.size() / 2);

    // open connectivity file:
    infile.open(conn_filename, std::ios::in);
    assert(infile.good() && "ReadMesh:: Connectivity file does not exist");
    connectivity.clear();

    // global node no of triangles:
    int node1, node2, node3;

    while (infile >> node1 >> node2 >> node3)
    {
        connectivity.push_back(node1);
        connectivity.push_back(node2);
        connectivity.push_back(node3);
    }

    infile.close();

    assert(connectivity.size() % 3 == 0);
    const int nElements = static_cast<int>(connectivity.size() / 3);

    std::cout << "\nRead " << nNodes << " nodes and " << nElements << " elements \n"
              << std::flush;
    return;
}

// write solution in tecplot format:
void write(const std::string filename,
                const std::vector<double> &coordinates,
                const std::vector<int> &connectivity,
                const double *nodal_field = nullptr)
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
    const int nNodes = static_cast<int>(coordinates.size() / 2);
    const int nElements = static_cast<int>(connectivity.size() / 3);
    outfile << "ZONE t=\"t:0\", N=" << nNodes << ", E=" << nElements << ", F=FEPOINT, ET=TRIANGLE\n";

    // Nodal coordinates and field values
    for (int i = 0; i < nNodes; ++i)
    {
        const double *X = &coordinates[2 * i];
        if (nodal_field != nullptr)
            outfile << X[0] << " " << X[1] << " " << nodal_field[i] << "\n";
        else
            outfile << X[0] << " " << X[1] << "\n";
    }

    // Connectivity
    for (int e = 0; e < nElements; ++e)
    {
        const int *conn = &connectivity[3 * e];
        outfile << conn[0] + 1 << " " << conn[1] + 1 << " " << conn[2] + 1 << "\n";
    }

    outfile.close();
}