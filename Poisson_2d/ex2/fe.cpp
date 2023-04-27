#include "fe.h"

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
        break;
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

