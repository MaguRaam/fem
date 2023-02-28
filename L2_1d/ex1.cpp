// compute L2 projection of a given function f -> Pf:
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cassert>
#include <functional>
#include <Eigen/Sparse>

std::array<std::function<double(double)>, 2> basis_function(double, double);
double simpsons_rule(std::function<double(double)>, double, double);
void local_force_vec(std::function<double(double)>, double, const double, double Flocal[2]);
void local_mass_matrix(const double, const double, double Mlocal[2][2]);
double local_L2_error(std::function<double(double)>, const double, const double, const double, const double);


int main()
{

    // define function f;
    auto f = [](double x)
    { return x; };

    // domain to define Vh space:
    const double L = 1.0;

    // no of Nodes:
    const int nNodes = 40;

    // no of Elements:
    const int nElements = nNodes - 1;

    // element size:
    const double dh = L / static_cast<double>(nElements);

    // define nodal coordinates:
    std::vector<double> x(nNodes);

    for (int n = 0; n < nNodes; ++n)
        x[n] = static_cast<double>(n) * dh;

    // connectivity matrix:
    std::vector<std::array<int, 2>> connectivity(nElements);

    for (int e = 0; e < nElements; ++e)
    {
        connectivity[e][0] = e;
        connectivity[e][1] = e + 1;
    }

    // Initialize sparse Mass matrix:
    Eigen::SparseMatrix<double, Eigen::RowMajor> M(nNodes, nNodes);
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nElements * 4);

    // Initalize force vector:
    Eigen::VectorXd F(nNodes);
    F.setZero();

    // Initilize solution vector:
    Eigen::VectorXd U(nNodes);
    U.setZero();

    // assemble force vector and mass matrix:
    for (int e = 0; e < nElements; ++e)
    {
        // get nodal coordinates of the element:
        double a = x[e];
        double b = x[e + 1];

        // initialize local force vector:
        double Flocal[2];

        // initialize local mass matrix:
        double Mlocal[2][2];

        // compute local force vector:
        local_force_vec(f, a, b, Flocal);

        // compute local mass matrix:
        local_mass_matrix(a, b, Mlocal);

        // add local force to global force vector:
        F(e) += Flocal[0];
        F(e + 1) += Flocal[1];

        // push back element mass matrix to tripletList:
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                tripletList.push_back({connectivity[e][i], connectivity[e][j], Mlocal[i][j]});
    }

    tripletList.shrink_to_fit();

    // Actual assembly of global mass matrix:
    M.setFromTriplets(tripletList.begin(), tripletList.end());

    // solve MU = F:
    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
    solver.analyzePattern(M);
    solver.factorize(M);
    U = solver.solve(F);
    assert(solver.info() == Eigen::Success);

    // write solution:
    std::ofstream file("sol.dat", std::ios::out);
    file.flags(std::ios::dec | std::ios::scientific);
    file.precision(16);

    for (int n = 0; n < nNodes; ++n)
        file << x[n] << "\t" << U[n] << "\t" << f(x[n]) << "\n";

    file.close();

    // compute L2 error:
    file.open("error.dat", std::ios::app);

    double L2_error = 0.0;
    
    for (int e = 0; e < nElements; ++e)
        L2_error += local_L2_error(f, x[e], x[e+1], U[e], U[e+1]);

    file << dh << "\t" << sqrt(L2_error) << "\n";

    return 0;
}

// compute local force vector:
void local_force_vec(std::function<double(double)> f, double a, double b, double Flocal[2])
{
    //basis:
    auto N = basis_function(a, b);

    // fill local force vector: Fi = (f, N[i]):
    for (int i = 0; i < 2; ++i){
    	
    	auto fNi = [&](double x){ return f(x)*N[i](x); };
    	Flocal[i] = simpsons_rule(fNi, a, b);

    }
}

// compute local mass matrix:
void local_mass_matrix(double a, double b, double Mlocal[2][2])
{
    // basis:
    auto N = basis_function(a, b);

    // fill local mass matrix: Mij = (N[i], N[j]):
    for (int i = 0; i < 2; ++i){
        for (int j = 0; j < 2; ++j){

            auto NiNj = [&](double x){return N[i](x) * N[j](x);};
            Mlocal[i][j] = simpsons_rule(NiNj, a, b);
        }
    }
    
}

// compute Local L2 error:
double local_L2_error(std::function<double(double)> f, const double a, const double b, const double U0, const double U1){

    // basis:
    auto N = basis_function(a, b);

    // construct Pf(x) = U0N[0](x) + U1N[1](x) on the element:
    auto Pf = [&](double x){return U0*N[0](x) + U1*N[1](x);};

    // error function e(x)^2 = (f(x) - Pf(x))^2;
    auto e_sqr = [&](double x){return pow(f(x) - Pf(x), 2.0);};

    // integrate square error:
    return simpsons_rule(e_sqr, a, b);
}


// basis function:
std::array<std::function<double(double)>, 2> basis_function(double a, double b){

    // basis function:
    std::array<std::function<double(double)>, 2> N;
    
    N[0] = [a, b](double x){ return (b - x) / (b - a);};
    N[1] = [a, b](double x){ return (x - a) / (b - a);};
    
    return N;
}


// integrate a function using simpson's rule:
double simpsons_rule(std::function<double(double)> f, double a, double b)
{
    return (b - a) * (f(a) + 4.0 * f(0.5 * (a + b)) + f(b)) * (1.0 / 6.0);
}

