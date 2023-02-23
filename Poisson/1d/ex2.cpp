// solve 1d poisson equation with nonzero dirclet boundary conditions:

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cassert>
#include <Eigen/Sparse>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

// compute local stiffness matrix Klocal:
inline Matrix2d k_local(const double &a, const double &b)
{
    double h = b - a;
    return Matrix2d{{1.0 / h, -1.0 / h}, {-1.0 / h, 1.0 / h}};
}

// compute local force vector Flocal:
inline Vector2d f_local(const double &a, const double &b)
{
    double h = b - a;
    return Vector2d{0.5 * h, 0.5 * h};
}

// integrate a function using simpson's rule:
inline double simpsons_rule(std::function<double(double)> f, const double &a, const double &b)
{
    return (b - a) * (f(a) + 4.0 * f(0.5 * (a + b)) + f(b)) * (1.0 / 6.0);
}

// compute local L2error:
inline double local_L2_error(std::function<double(double)> f, const double &a, const double &b, const double &Ua, const double &Ub)
{

    auto Na = [&](double x)
    { return (b - x) / (b - a); };
    auto Nb = [&](double x)
    { return (x - a) / (b - a); };
    auto e_sqr = [&](double x)
    { return pow(f(x) - Ua * Na(x) - Ub * Nb(x), 2.0); };

    return simpsons_rule(e_sqr, a, b);
}

int main()
{
    // create grid points to define Vh space on [0, 1] domain:
    const int nNodes = 2000;
    const int nElements = nNodes - 1;
    const double dh = 1.0 / static_cast<double>(nElements);

    vector<double> x(nNodes);
    for (int n = 0; n < nNodes; ++n)
        x[n] = n * dh;

    // compute exact solution:
    VectorXd Uexact(nNodes);
    auto uexact = [](double xi)
    { return 0.5 * xi * (1 - xi) + 1.0; };
    std::transform(x.begin(), x.end(), Uexact.begin(), uexact);

    // connectivity matrix:
    MatrixX2i conn(nElements, 2);
    for (int e = 0; e < nElements; ++e)
    {
        conn(e, 0) = e;
        conn(e, 1) = e + 1;
    }

    // initialize stiffness matrix and force vector:
    SparseMatrix<double, Eigen::RowMajor> K(nNodes, nNodes);
    VectorXd F = VectorXd::Zero(nNodes);

    // Assemble Sitffness matrix K and Force vector F for only interior nodes:
    vector<Triplet<double>> tripleList;
    tripleList.reserve(nElements * 4);

    // loop over interior elements:
    for (int e = 0; e < nElements; ++e)
    {
        // get nodal values:
        const double &a = x[e];
        const double &b = x[e + 1];

        // get local stiffness matrix and force vector:
        auto kLocal = k_local(a, b);
        auto fLocal = f_local(a, b);

        // loop over nodes in the element:
        for (int i = 0; i < 2; ++i)
        {
            const int &I = conn(e, i); // get global node number:

            for (int j = 0; j < 2; ++j)
            {
                const int &J = conn(e, j); // get global node number:

                tripleList.push_back({I, J, kLocal(i, j)});
            }

            // force vector:
            F(I) += fLocal(i);
        }
    }

    // build sparse matrix K:
    K.setFromTriplets(tripleList.begin(), tripleList.end());

    //Enforce boundary conditions:

    // make first and last rows zero:
    for (int j = 0; j < nNodes; ++j){
        K.coeffRef(0, j) = 0.0;
        K.coeffRef(nNodes - 1, j) = 0.0;
    }

    // set diagonal elements as 1:
    K.coeffRef(0, 0) = K.coeffRef(nNodes - 1, nNodes - 1) = 1.0;

    //! set boundary condition by modifying force vector:
    F.coeffRef(0) = 1.0;
    F.coeffRef(nNodes - 1) = 1.0;


    // solve MU = F:
    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
    solver.analyzePattern(K);
    solver.factorize(K);
    auto U = solver.solve(F);
    assert(solver.info() == Eigen::Success);

    // write solution:
    std::ofstream file("sol.dat", std::ios::out);
    file.flags(std::ios::dec | std::ios::scientific);
    file.precision(16);

    for (int n = 0; n < nNodes; ++n)
        file << x[n] << " " << U(n) << " " << Uexact(n) << "\n";

    file.close();

     // compute L2 error:
    file.open("error.dat", std::ios::app);

    // square of error function e(x) = u(x) - uh(x):
    double L2_error = 0.0;
    for (int e = 1; e < nElements - 1; ++e)
        L2_error += local_L2_error(uexact, x[e], x[e + 1], U(conn(e, 0)), U(conn(e, 1))); 

    file << dh << "\t" << sqrt(L2_error) << "\n";

}