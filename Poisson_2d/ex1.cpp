// Solve Poisson equation with Dirchlet boundary condition on a triangular mesh:

#include <Eigen/Sparse>
#include <iostream>
#include "io.h"
#include "fe.h"
#include "tools.h"

using namespace std;
using namespace Eigen;

// compute element stiffness matrix:
inline Matrix3d compute_Klocal(const Vector3d &X, const Vector3d &Y)
{
    Matrix3d Klocal = Matrix3d::Zero();

    // get quadrature pts and weights:
    const auto &nqpts = GaussQuadratureTriangle<1>::nqpts;
    const auto &wts = GaussQuadratureTriangle<1>::wts;
    const auto &r = GaussQuadratureTriangle<1>::r;
    const auto &s = GaussQuadratureTriangle<1>::s;

    // loop over quadrature pts and compute local stiffness matrix:
    for (int q = 0; q < nqpts; ++q)
    {
        // get basis gradients and Jacobian at quadpoint:
        auto [S, dSdx, dSdy, detJ] = isomap(r[q], s[q], X, Y);

        Klocal += (dSdx * dSdx.transpose() + dSdy * dSdy.transpose()) * detJ * 0.5 * wts[q]; //! note we scale weight by 0.5:
    }

    return Klocal;
}

// compute element force vector:
template <typename F>
inline Vector3d compute_Flocal(F f, const Vector3d &X, const Vector3d &Y)
{
    Vector3d Flocal = Vector3d::Zero();

    // get quadrature pts and weights:
    const auto &nqpts = GaussQuadratureTriangle<1>::nqpts;
    const auto &wts = GaussQuadratureTriangle<1>::wts;
    const auto &r = GaussQuadratureTriangle<1>::r;
    const auto &s = GaussQuadratureTriangle<1>::s;

    // loop over quadrature pts and compute local force vector:
    for (int q = 0; q < nqpts; ++q)
    {
        // get basis, basis gradients and Jacobian at quadpoint:
        auto [S, dSdx, dSdy, detJ] = isomap(r[q], s[q], X, Y);

        // map r, s -> x, y:
        double x = X.dot(S);
        double y = Y.dot(S);

        Flocal += f(x, y) * S * detJ * 0.5 * wts[q];
    }

    return Flocal;
}

int main()
{
    // exact solution:
    auto uexact = [](const double &x, const double &y)
    { return sin(M_PI * x) * sin(M_PI * y); };

    // gradient of exact solution:
    auto grad_uexact = [](const double &x, const double &y)
    { return Vector2d(M_PI * cos(M_PI * x) * sin(M_PI * y), M_PI * cos(M_PI * y) * sin(M_PI * x)); };

    // forcing function:
    auto f = [](const double &x, const double &y)
    { return 2.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y); };

    // load connectivity matrix:
    MatrixX3i conn;
    load_connectivity("mesh/ele.dat", conn);

    // load nodal coordinates:
    VectorXd x, y;
    load_coordinates("mesh/node.dat", x, y);

    // no of elements and nodes:
    const int nElements = conn.rows();
    const int nNodes = x.size();

    // compute exact solution:
    VectorXd Uexact(nNodes);
    for (int n = 0; n < nNodes; ++n)
        Uexact(n) = uexact(x(n), y(n));

    // assemble stiffness matrix K and force vector F:
    SparseMatrix<double, RowMajor> K(nNodes, nNodes);
    VectorXd F = VectorXd::Zero(nNodes);

    // create tripletList:
    vector<Triplet<double>> tripleList;
    tripleList.reserve(nElements * 9);

    // loop over elements:
    for (int e = 0; e < nElements; ++e)
    {
        // get coordinates of the triangle nodes:
        Vector3d X, Y;
        for (int n = 0; n < 3; ++n)
        {
            X(n) = x(conn(e, n));
            Y(n) = y(conn(e, n));
        }

        // get local stiffness matrix and force vector given vertex coordinates of the element:
        Matrix3d Klocal = compute_Klocal(X, Y);
        Vector3d Flocal = compute_Flocal(f, X, Y);

        // loop over nodes in the element:
        for (int i = 0; i < 3; ++i)
        {
            // push back matrix:
            for (int j = 0; j < 3; ++j)
                tripleList.emplace_back(conn(e, i), conn(e, j), Klocal(i, j));

            // assemble force:
            F(conn(e, i)) += Flocal(i);
        }
    }

    // build sparse matrix K from triplet list:
    K.setFromTriplets(tripleList.begin(), tripleList.end());

    // set dirchlet boundary condition:
    for (int n = 0; n < nNodes; ++n)
    {
        // check if node is located on the boundary:
        if (x(n) == 0.0 || x(n) == 1.0 || y(n) == 0.0 || y(n) == 1.0)
        {
            // make the nth row zeros:
            for (int j = 0; j < nNodes; ++j)
                K.coeffRef(n, j) = 0.0;

            // make the diagonal element 1:
            K.coeffRef(n, n) = 1.0;

            // set boundary conditon
            F(n) = uexact(x(n), y(n));
        }
    }

    // solve KU = F:
    SparseLU<decltype(K)> solver;
    solver.analyzePattern(K);
    solver.factorize(K);
    VectorXd U = solver.solve(F);
    assert(solver.info() == Success);

    // plot:

    // compute abs Error:
    VectorXd Error = (U - Uexact).cwiseAbs();

    // write solution in .vtk:
    write_vtk(conn, x, y, U, Uexact, Error, "sol.vtk");

    // Postprocessing:

    // compute l2error, energy error, energy functional of uexact, uh and piu:
    double l2error = 0.0;
    double energyerror_uh = 0.0;
    double energyerror_piu = 0.0;
    double energy_uexact = 0.0;
    double energy_uh = 0.0;
    double energy_piu = 0.0;

    for (int e = 0; e < nElements; ++e)
    {
        // get coordinates of the nodes and the nodal values:
        Vector3d X, Y, Uh, piU;
        for (int n = 0; n < 3; ++n)
        {
            X(n) = x(conn(e, n));
            Y(n) = y(conn(e, n));
            Uh(n) = U(conn(e, n));
            piU(n) = Uexact(conn(e, n));
        }

        energy_uh += energy_functional(X, Y, Uh, f);
        energy_piu += energy_functional(X, Y, piU, f);
        energy_uexact += energy_functional(X, Y, uexact, grad_uexact, f);

        l2error += L2error(uexact, X, Y, Uh);
        energyerror_uh += energy_error(grad_uexact, X, Y, Uh);
        energyerror_piu += energy_error(grad_uexact, X, Y, piU);
    }

    // write L2 and Energy error to a file:
    std::ofstream file("error.dat", std::ios::app);
    file.flags(std::ios::dec | std::ios::scientific);

    file << sqrt(1.0 / nElements) << "\t" << sqrt(l2error) << "\t" << sqrt(energyerror_uh) << std::endl;
    file.close();

    // write energy of exact, fem and interpolated solution:
    file.open("energy.dat", std::ios::app);
    file << sqrt(1.0 / nElements) << "\t" << energy_uexact << "\t" << energy_uh << "\t" << energy_piu << "\n";
    file.close();

    // write energy error of fem and interpolated solution:
    file.open("energy_error.dat", std::ios::app);
    file << sqrt(1.0 / nElements) << "\t" << energyerror_uh << "\t" << energyerror_piu << "\n";
    file.close();

    return 0;
}
