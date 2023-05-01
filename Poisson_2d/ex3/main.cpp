#include <iostream>
#include <vector>
#include <cmath>

#include "fe.h"

// forcing function:
double f(double x, double y) { return 2.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y); }

// exact solution:
double uexact(double x, double y) { return sin(M_PI * x) * sin(M_PI * y); }

// element stiffness matrix:
void GetElementStiffnessMatrix(const std::vector<double> &nodal_coords, double kmat[4][4]);

// element force vector:
void GetElementForceVector(const std::vector<double> &nodal_coords, double *fvec, double (*f)(double, double));

int main()
{
    // map parametric element to physical element.
    std::vector<double> nodal_coords{-1., -1., 1., -1., 1., 1., -1., 1.};

    return 0;
}

// element stiffness matrix:
void GetElementStiffnessMatrix(const std::vector<double> &nodal_coords, double kmat[4][4])
{
    // loop over quadrature points:
    for (int q = 0; q < nqpts; ++q)
    {
        double detJ;

        // fill matrix:
        for (int i = 0; i < 4; ++i)
        {
            double gradNi[2], Ni;
            GetShapeFunction(xi[q], i, nodal_coords, Ni, gradNi, detJ);

            for (int j = 0; j < 4; ++j)
            {
                double gradNj[2], Nj;
                GetShapeFunction(xi[q], j, nodal_coords, Nj, gradNj, detJ);

                kmat[i][j] += (gradNi[0] * gradNj[0] + gradNi[1] * gradNj[1]) * detJ * wts[q];
            }
        }
    }
}

// element force vector:
void GetElementForceVector(const std::vector<double> &nodal_coords, double *fvec, double (*f)(double, double))
{
    // loop over quadrature points:
    for (int q = 0; q < nqpts; ++q)
    {
        // map point in parametric element to physical element:
        double detJ, x[2], J[2][2];
        ComputeIsoparamMap(xi[q], nodal_coords, x, J);

        // fill force vector:
        for (int i = 0; i < 4; ++i)
        {
            double Ni, gradNi[2];
            GetShapeFunction(xi[q], i, nodal_coords, Ni, gradNi, detJ);
            fvec[i] += f(x[0], x[1]) * Ni * detJ * wts[q];  
        }
    }
}