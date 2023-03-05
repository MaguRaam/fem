#include <Eigen/Sparse>
#include "io.h"

using namespace std;
using namespace Eigen;

// quadarture pts and weights:
static const int nqpts = 1;
static const double wts[nqpts] = {1.0};
static const double xi[nqpts][2] = {{1.0 / 3.0, 1.0 / 3.0}};

// shape functions and derivatives over the parametric element:
void GetShapeFunction(const double *xi, int a, double &val, double *dval);

// shape functions and derivatives over the physical element:
void GetShapeFunction(const double *xi, int a, const double *nodal_coords, double &val, double *dval, double &detJ);

// element stiffness matrix:
void GetElementStiffnessMatrix(const double *nodal_coords, double kmat[][3]);

// element force vector:
void GetElementForceVector(const double *nodal_coords, double *fvec, double (*f)(double, double));

// compute L2 error:
double L2Error(const double *nodal_coords, const double *unode, double (*uexact)(double, double));

// compute energy error:
double EnergyError(const double *nodal_coords, const double *unode, void (*grad_uexact)(double, double, double *graduexact));

// compute energy functional of fem or interpolated solution:
double EnergyFunctional(const double *nodal_coords, const double *unode, double (*f)(double, double));

// compute energy functional of fem or interpolated solution:
double EnergyFunctional(const double *nodal_coords, double (*uexact)(double, double), void (*grad_uexact)(double, double, double *graduexact), double (*f)(double, double));

// forcing function:
double f(double x, double y) { return 2.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y); }

// exact solution:
double uexact(double x, double y) { return sin(M_PI * x) * sin(M_PI * y); }

// gradien of exact solution:
void grad_uexact(double x, double y, double graduexact[2])
{
    graduexact[0] = M_PI * cos(M_PI * x) * sin(M_PI * y);
    graduexact[1] = M_PI * cos(M_PI * y) * sin(M_PI * x);
}

int main()
{
    // read coordinates and connectivity from mesh file:
    std::vector<double> coordinates;
    std::vector<int> connectivity;
    read_mesh("mesh/node.dat", coordinates, "mesh/ele.dat", connectivity);

    // no of nodes and elements:
    const int nNodes = static_cast<int>(coordinates.size() / 2);
    const int nElements = static_cast<int>(connectivity.size() / 3);

    // assemble stiffness matrix K and force vector F:
    SparseMatrix<double, RowMajor> K(nNodes, nNodes);
    VectorXd F = VectorXd::Zero(nNodes);

    // create tripletList:
    vector<Triplet<double>> tripleList;
    tripleList.reserve(nElements * 9);

    // loop over elements:
    for (int e = 0; e < nElements; ++e)
    {
        // get connectivity:
        const int *conn = &connectivity[3 * e];

        // get nodal coordinates:
        double nodal_coord[6];
        for (int n = 0; n < 3; ++n)
        {
            nodal_coord[2 * n] = coordinates[2 * conn[n]];         // x
            nodal_coord[2 * n + 1] = coordinates[2 * conn[n] + 1]; // y
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
    for (int n = 0; n < nNodes; ++n)
    {
        const double *x = &coordinates[2 * n];

        // check if node is located on the boundary:
        if (x[0] == 0.0 || x[0] == 1.0 || x[1] == 0.0 || x[1] == 1.0)
        {
            // make the nth row zeros:
            for (int j = 0; j < nNodes; ++j)
                K.coeffRef(n, j) = 0.0;

            // make the diagonal element 1:
            K.coeffRef(n, n) = 1.0;

            // set boundary conditon
            F(n) = uexact(x[0], x[1]);
        }
    }

    // solve KU = F:
    SparseLU<decltype(K)> solver;
    solver.analyzePattern(K);
    solver.factorize(K);
    VectorXd U = solver.solve(F);
    assert(solver.info() == Success);

    // compute exact solution:
    VectorXd Uexact(nNodes);
    for (int n = 0; n < nNodes; ++n)
    {
        const double *x = &coordinates[2 * n];
        Uexact(n) = uexact(x[0], x[1]);
    }

    // compute abs Error:
    VectorXd Error = (U - Uexact).cwiseAbs();

    // write solution:
    write("sol.dat", coordinates, connectivity, Error.data());

    // Postprocessing:

    // compute l2error:
    double l2error = 0.0, energyerror_uh = 0.0, energyerror_piu = 0.0;
    double energy_uexact = 0.0, energy_uh = 0.0, energy_piu = 0.0;

    for (int e = 0; e < nElements; ++e)
    {
        // get connectivity:
        const int *conn = &connectivity[3 * e];

        // get nodal coordinates:
        double nodal_coord[6], Uh[3], piU[3];
        for (int n = 0; n < 3; ++n)
        {
            nodal_coord[2 * n] = coordinates[2 * conn[n]];         // x
            nodal_coord[2 * n + 1] = coordinates[2 * conn[n] + 1]; // y
            Uh[n] = U(conn[n]);
            piU[n] = Uexact(conn[n]);
        }

        l2error += L2Error(nodal_coord, Uh, uexact);
        energyerror_uh += EnergyError(nodal_coord, Uh, grad_uexact);
        energyerror_piu += EnergyError(nodal_coord, piU, grad_uexact);
    
        energy_uh += EnergyFunctional(nodal_coord, Uh, f);
        energy_piu += EnergyFunctional(nodal_coord, piU, f);
        energy_uexact += EnergyFunctional(nodal_coord, uexact, grad_uexact, f);
    }

    // write L2 and Energy error to a file:
    std::ofstream file("error.dat", std::ios::app);
    file.flags(std::ios::dec | std::ios::scientific);
    file.precision(16);

    file << sqrt(1.0 / nElements) << "\t" << sqrt(l2error) << "\t" << sqrt(energyerror_uh) << std::endl;
    file.close();

    // write energy error of fem and interpolated solution:
    file.open("energy_error.dat", std::ios::app);
    file << sqrt(1.0 / nElements) << "\t" << energyerror_uh << "\t" << energyerror_piu << "\n";
    file.close();

    // write energy of exact, fem and interpolated solution:
    file.open("energy.dat", std::ios::app);
    file << sqrt(1.0 / nElements) << "\t" << energy_uexact << "\t" << energy_uh << "\t" << energy_piu << "\n";
    file.close();

    return 0;
}

// get shape and derivatives over the parametric element:
void GetShapeFunction(const double *xi, int a, double &val, double *dval)
{
    switch (a)
    {
    case 0:
        val = 1 - xi[0] - xi[1];
        dval[0] = -1.0;
        dval[1] = -1.0;
        break;
    case 1:
        val = xi[0];
        dval[0] = 1.0;
        dval[1] = 0.0;
        break;
    case 2:
        val = xi[1];
        dval[0] = 0.0;
        dval[1] = 1.0;
    default:
        break;
    }
}

// map parametric element to physical element:
void isomap(const double *xi, double *x, const double *nodal_coords)
{
    x[0] = 0.0;
    x[1] = 0.0;
    double Ni, gradNi[2];

    for (int i = 0; i < 3; ++i)
    {
        GetShapeFunction(xi, i, Ni, gradNi);

        x[0] += nodal_coords[2 * i] * Ni;
        x[1] += nodal_coords[2 * i + 1] * Ni;
    }
}

// shape functions and derivatives over the physical element:
void GetShapeFunction(const double *xi, int a, const double *nodal_coords, double &val, double *dval, double &detJ)
{
    // get gradient of shape function wrt xi:
    double dval_xi[2];
    GetShapeFunction(xi, a, val, dval_xi);

    // get vertex coordinates:
    double x0 = nodal_coords[0], y0 = nodal_coords[1];
    double x1 = nodal_coords[2], y1 = nodal_coords[3];
    double x2 = nodal_coords[4], y2 = nodal_coords[5];

    // compute transpose of Jacobian matrix:
    double JT[2][2] = {{x1 - x0, y1 - y0},
                       {x2 - x0, y2 - y0}};

    // compute detJ:
    detJ = JT[0][0] * JT[1][1] - JT[0][1] * JT[1][0];

    // inverse of the transpose of Jacobian matrix:
    double JTinv[2][2];

    JTinv[0][0] = JT[1][1] / detJ;
    JTinv[1][1] = JT[0][0] / detJ;
    JTinv[0][1] = -JT[0][1] / detJ;
    JTinv[1][0] = -JT[1][0] / detJ;

    // compute gradient wrt x,y:
    dval[0] = JTinv[0][0] * dval_xi[0] + JTinv[0][1] * dval_xi[1];
    dval[1] = JTinv[1][0] * dval_xi[0] + JTinv[1][1] * dval_xi[1];
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

// compute L2 error:
double L2Error(const double *nodal_coords, const double *unode, double (*uexact)(double, double))
{
    double l2error = 0.0;

    // loop over quadrature points:
    for (int q = 0; q < nqpts; ++q)
    {
        // map point in parametric element to physical element:
        double detJ, x[2];
        isomap(xi[q], x, nodal_coords);

        double uh = 0.0;

        // get shape function and detJ at xi:
        for (int i = 0; i < 3; ++i)
        {
            double Ni, gradNi[2];
            GetShapeFunction(xi[q], i, nodal_coords, Ni, gradNi, detJ);
            uh += Ni * unode[i];
        }

        l2error += pow((uexact(x[0], x[1]) - uh), 2.0) * detJ * 0.5 * wts[q];
    }

    return l2error;
}

// compute energy error:
double EnergyError(const double *nodal_coords, const double *unode, void (*grad_uexact)(double, double, double *graduexact))
{
    double error = 0.0;

    // loop over quadrature points:
    for (int q = 0; q < nqpts; ++q)
    {
        // map point in parametric element to physical element:
        double detJ, x[2];
        isomap(xi[q], x, nodal_coords);

        double grad_uh[2] = {0.0, 0.0};

        // get shape function and detJ at xi and compute gradient of uh:
        for (int i = 0; i < 3; ++i)
        {
            double Ni, gradNi[2];
            GetShapeFunction(xi[q], i, nodal_coords, Ni, gradNi, detJ);
            grad_uh[0] += gradNi[0] * unode[i];
            grad_uh[1] += gradNi[1] * unode[i];
        }

        // get exact gradient:
        double graduexact[2];
        grad_uexact(x[0], x[1], graduexact);

        error += ((graduexact[0] * graduexact[0] + graduexact[1] * graduexact[1]) + (grad_uh[0] * grad_uh[0] + grad_uh[1] * grad_uh[1]) - 2.0 * (graduexact[0] * grad_uh[0] + graduexact[1] * grad_uh[1])) * detJ * 0.5 * wts[q];
    }

    return error;
}

// compute energy functional of fem or interpolated solution:
double EnergyFunctional(const double *nodal_coords, const double *unode, double (*f)(double, double))
{
    double energy = 0.0;

    // loop over quadrature points:
    for (int q = 0; q < nqpts; ++q)
    {
        // map point in parametric element to physical element:
        double detJ, x[2];
        isomap(xi[q], x, nodal_coords);

        double grad_u[2] = {0.0, 0.0}, u = 0.0;

        // get shape function and detJ at xi and compute gradient of uh:
        for (int i = 0; i < 3; ++i)
        {
            double Ni, gradNi[2];
            GetShapeFunction(xi[q], i, nodal_coords, Ni, gradNi, detJ);
            grad_u[0] += gradNi[0] * unode[i];
            grad_u[1] += gradNi[1] * unode[i];
            u += Ni * unode[i];
        }

        energy += (0.5 * (grad_u[0] * grad_u[0] + grad_u[1] * grad_u[1]) - f(x[0], x[1]) * u) * detJ * 0.5 * wts[q];
    }

    return energy;
}

// compute energy functional of fem or interpolated solution:
double EnergyFunctional(const double *nodal_coords, double (*uexact)(double, double), void (*grad_uexact)(double, double, double *graduexact), double (*f)(double, double))
{
    double energy = 0.0;

    // loop over quadrature points:
    for (int q = 0; q < nqpts; ++q)
    {
        // map point in parametric element to physical element:
        double detJ, x[2];
        isomap(xi[q], x, nodal_coords);

        // compute detJ
        double Ni, gradNi[2];
        GetShapeFunction(xi[q], 0, nodal_coords, Ni, gradNi, detJ);

        // get exact gradient:
        double grad_u[2];
        grad_uexact(x[0], x[1], grad_u);
        double u = uexact(x[0], x[1]);
    
        energy += (0.5 * (grad_u[0] * grad_u[0] + grad_u[1] * grad_u[1]) - f(x[0], x[1])*u) * detJ * 0.5 * wts[q];
    }

    return energy;
}