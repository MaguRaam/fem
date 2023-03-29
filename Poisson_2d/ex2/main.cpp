#include <Eigen/Sparse>
#include <iostream>
#include <fstream>
#include <cmath>

#include "mesh.h"
#include "fe.h"

using namespace std;
using namespace Eigen;

// forcing function:
double f(double x, double y) { return 2.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y); }

// exact solution:
double uexact(double x, double y) { return sin(M_PI * x) * sin(M_PI * y); }

// element stiffness matrix:
void GetElementStiffnessMatrix(const double *nodal_coords, double kmat[][3]);

// element force vector:
void GetElementForceVector(const double *nodal_coords, double *fvec, double (*f)(double, double));

// set boundary conditon:
void SetDirichletBCs(const Mesh &M, SparseMatrix<double, RowMajor> &K, VectorXd &F);

int main()
{
    Mesh M;
    ReadMesh("mesh/coordinates.dat", "mesh/connectivity.dat", "mesh/boundarynodes.dat", M);

    // assemble stiffness matrix K and force vector F:
    SparseMatrix<double, RowMajor> K(M.nNodes, M.nNodes);
    VectorXd F = VectorXd::Zero(M.nNodes);

    // create tripletList:
    vector<Triplet<double>> tripleList;
    tripleList.reserve(M.nElements * 9);

    // loop over elements:
    for (int e = 0; e < M.nElements; ++e)
    {
        // get connectivity:
        const int *conn = &M.connectivity[3 * e];

        // get nodal coordinates:
        double nodal_coord[6];
        for (int n = 0; n < 3; ++n)
        {
            nodal_coord[2 * n] = M.coordinates[2 * conn[n]];         // x
            nodal_coord[2 * n + 1] = M.coordinates[2 * conn[n] + 1]; // y
        }

        // initalize local element matrix:
        double kmat[3][3] = {{0.0}};
        double fvec[3] = {0.0};

        GetElementStiffnessMatrix(nodal_coord, kmat);
        GetElementForceVector(nodal_coord, fvec, f);

        // add local stiffness matrix and force vector to global:
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                tripleList.emplace_back(conn[i], conn[j], kmat[i][j]);
            }
            F.coeffRef(conn[i]) += fvec[i];
        }
    }

    // build sparse matrix K from triplet list:
    K.setFromTriplets(tripleList.begin(), tripleList.end());

    // set dirchlet boundary condition:
    SetDirichletBCs(M, K, F);

    // solve KU = F:
    SparseLU<decltype(K)> solver;
    solver.analyzePattern(K);
    solver.factorize(K);
    VectorXd U = solver.solve(F);
    assert(solver.info() == Success);

    WriteSolution("sol.dat", M, U.data());

    // compute exact solution:
    VectorXd Uexact(M.nNodes);
    for (int n = 0; n < M.nNodes; ++n)
    {
        const double *x = &M.coordinates[2 * n];
        Uexact(n) = uexact(x[0], x[1]);
    }

    // write L2 and Energy error to a file:
    std::ofstream file("error.dat", std::ios::app);
    file.flags(std::ios::dec | std::ios::scientific);
    file.precision(16);

    double l2error = (1.0/sqrt(M.nNodes))*(U - Uexact).norm();

    file << sqrt(1.0 / M.nElements) << "\t" << l2error << "\n";

    file.close();

    return 0;
}

// element stiffness matrix:
void GetElementStiffnessMatrix(const double *nodal_coords, double kmat[][3])
{

    // loop over quadrature points:
    for (int q = 0; q < nqpts; ++q)
    {
        double detJ;

        // fill matrix:
        for (int i = 0; i < 3; ++i)
        {
            double gradNi[2], Ni;
            GetShapeFunction(xi[q], i, nodal_coords, Ni, gradNi, detJ);

            for (int j = 0; j < 3; ++j)
            {
                double gradNj[2], Nj;
                GetShapeFunction(xi[q], j, nodal_coords, Nj, gradNj, detJ);

                kmat[i][j] += (gradNi[0] * gradNj[0] + gradNi[1] * gradNj[1]) * 0.5 * detJ * wts[q];
            }
        }
    }
}

// element force vector:
void GetElementForceVector(const double *nodal_coords, double *fvec, double (*f)(double, double))
{
    // loop over quadrature points:
    for (int q = 0; q < nqpts; ++q)
    {
        // map point in parametric element to physical element:
        double detJ, x[2];
        isomap(xi[q], x, nodal_coords);

        // fill force vector:
        for (int i = 0; i < 3; ++i)
        {
            double Ni, gradNi[2];
            GetShapeFunction(xi[q], i, nodal_coords, Ni, gradNi, detJ);
            fvec[i] += f(x[0], x[1]) * Ni * 0.5 * detJ * wts[q];
        }
    }
}

void SetDirichletBCs(const Mesh &M, SparseMatrix<double, RowMajor> &K, VectorXd &F)
{
    for (const auto &n : M.bdNodes)
    {
        // make the nth row zeros:
        for (int j = 0; j < M.nNodes; ++j)
            K.coeffRef(n, j) = 0.0;
        
        // make the diagonal element 1:
        K.coeffRef(n, n) = 1.0;

        // set boundary conditon
        F(n) = 0.0;
    }
}