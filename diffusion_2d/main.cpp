#include "mesh.h"
#include <iostream>
#include <cmath>

int main()
{
    Mesh M;
    ReadMesh("mesh/coordinates.dat", "mesh/connectivity.dat", "mesh/boundarynodes.dat", M);

    std::vector<double> U(M.nNodes);
    for (int n = 0; n < M.nNodes; ++n){
        double x = M.coordinates[2*n];
        double y = M.coordinates[2*n + 1];
        U[n] = sin(2*M_PI*x)*sin(2*M_PI*y);        
    }

    WriteSolution("sol.dat", M, U.data());

    return 0;
}
