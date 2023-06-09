#pragma once

#include <Eigen/Dense>
#include <fstream>

// load connectivity matrix for triangular mesh:
void load_connectivity(const std::string& filename, Eigen::MatrixX3i &conn)
{
    std::ifstream file(filename);
    if (!file.good())
        throw std::runtime_error("Failed to open " + filename);

    int nElements;
    file >> nElements;

    conn.resize(nElements, 3);

    for (int e = 0; e < nElements; ++e)
        file >> conn(e, 0) >> conn(e, 1) >> conn(e, 2);
    file.close();
}

// load nodal coordinates:
void load_coordinates(const std::string& filename, Eigen::VectorXd &x, Eigen::VectorXd &y)
{
    std::ifstream file(filename);
    if (!file.good())
        throw std::runtime_error("Failed to open " + filename);

    int nNodes;
    file >> nNodes;

    x.resize(nNodes);
    y.resize(nNodes);

    for (int n = 0; n < nNodes; ++n)
        file >> x(n) >> y(n);
    file.close();
}

// write solution in vtk file format:
void write_vtk(const Eigen::MatrixX3i& conn, const Eigen::VectorXd& x, const Eigen::VectorXd& y, const Eigen::VectorXd& U, const Eigen::VectorXd& Uexact, const Eigen::VectorXd& Error, const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.good())
        throw std::runtime_error("Failed to open file " + filename);

    file.flags(std::ios::dec | std::ios::scientific);
    file.precision(16);

    // Write VTK header
    file << "# vtk DataFile Version 3.0\nvtk fileput\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    
    
    const int nNodes = x.size();
    file << "POINTS " << nNodes << " float\n";
    for (int n = 0; n < nNodes; ++n)
        file << x(n) << " " << y(n) << " 0\n";
    
    const int nElements = conn.rows();
    file << "CELLS " << nElements << " " << nElements * (3 + 1) << "\n";
    for (int e = 0; e < nElements; ++e)
        file << "3 " << conn(e, 0) << " " << conn(e, 1) << " " << conn(e, 2) << "\n";
    
    file << "CELL_TYPES " << nElements << "\n";
    for (int e = 0; e < nElements; ++e)
        file << "5\n";
    
    file << "POINT_DATA " << nNodes << "\n";
    
    file << "SCALARS U float 1\nLOOKUP_TABLE default\n";
    for (int n = 0; n < nNodes; ++n)
        file << U(n) << "\n";
    
    file << "SCALARS Uexact float 1\nLOOKUP_TABLE default\n";
    for (int n = 0; n < nNodes; ++n)
        file << Uexact(n) << "\n";
    
    file << "SCALARS Error float 1\nLOOKUP_TABLE default\n";
    for (int n = 0; n < nNodes; ++n)
        file << Error(n) << "\n";


    file.close();
}

