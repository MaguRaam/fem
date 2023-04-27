#include "fe.h"

void GetShapeFunction(const double *xi, // in:  parametric coordinates over Khat
                      const int a,      // in:  basis function number = 0,1,2 or 3
                      double &Na,       // out: basis function value Na
                      double *dNa)      // out: gradient of the basis function grad_Na
{
    assert(a >= 0 && a <= 3 && "Invalid value of 'a': must be 0, 1, 2, or 3");

    switch (a)
    {
    case 0:
        Na = 0.25 * (1 - xi[0]) * (1 - xi[1]);
        dNa[0] = -0.25 * (1 - xi[1]);
        dNa[1] = -0.25 * (1 - xi[0]);
        break;
    case 1:
        Na = 0.25 * (1 + xi[0]) * (1 - xi[1]);
        dNa[0] = 0.25 * (1 - xi[1]);
        dNa[1] = -0.25 * (1 + xi[0]);
        break;
    case 2:
        Na = 0.25 * (1 + xi[0]) * (1 + xi[1]);
        dNa[0] = 0.25 * (1 + xi[1]);
        dNa[1] = 0.25 * (1 + xi[0]);
        break;
    case 3:
        Na = 0.25 * (1 - xi[0]) * (1 + xi[1]);
        dNa[0] = -0.25 * (1 + xi[1]);
        dNa[1] = 0.25 * (1 - xi[0]);
        break;

    default:
        break;
    }
}

// compute iso parametric map and Jacobian of the map:
void ComputeIsoparamMap(const double *xi,                        // in: parametric coordinates in Khat
                        const std::vector<double> &nodal_coords, // in: nodal coordinates in K
                        double phi[2],                           // out: isoparametric map
                        double grad_phi[2][2])                   // out: gradient of the isoparametric map
{
    // compute x, y for given psi and eta:
    phi[0] = 0.0;
    phi[1] = 0.0;
    grad_phi[0][0] = grad_phi[0][1] = grad_phi[1][0] = grad_phi[1][1] = 0.0;

    // x = N0x0 + ... + N3x3, y = N0y0 + ... + N3y3
    for (int a = 0; a < 4; ++a)
    {
        double Na;
        double dNa[2];
        GetShapeFunction(xi, a, Na, dNa);

        // compute isoparametric map:
        phi[0] += nodal_coords[2 * a] * Na;
        phi[1] += nodal_coords[2 * a + 1] * Na;

        // fill jacobian matrix:
        grad_phi[0][0] += nodal_coords[2 * a] * dNa[0];
        grad_phi[0][1] += nodal_coords[2 * a] * dNa[1];
        grad_phi[1][0] += nodal_coords[2 * a + 1] * dNa[0];
        grad_phi[1][1] += nodal_coords[2 * a + 1] * dNa[1];
    }
}

// determinant of Jacobian:
double Determinant(const double grad_phi[2][2])
{
    return grad_phi[0][0] * grad_phi[1][1] - grad_phi[0][1] * grad_phi[1][0];
}

// inverse of Jacobian transpose:
void InverseJacTranspose(const double J[2][2], // in:  Jacobian of map
                         double JTinv[2][2]    // out: inverse of Jacobian transpose
)
{
    double detJ = Determinant(J);

    JTinv[0][0] = J[1][1] / detJ;
    JTinv[0][1] = -J[1][0] / detJ;
    JTinv[1][0] = -J[0][1] / detJ;
    JTinv[1][1] = J[0][0] / detJ;
}

// shape functions and derivatives over the physical element:
void GetShapeFunction(const double *xi,                        // in:  parametric coordinates over Khat
                      int a,                                   // in:  basis function number = 0,1,2 or 3
                      const std::vector<double> &nodal_coords, // in: nodal coordinates in K
                      double &val,                             // out: basis function value
                      double *dval,                            // out: gradient of the basis function grad_Na
                      double &detJ                             // determinant of Jacobian
)
{
    // get gradient of shape function wrt xi:
    double dval_xi[2];
    GetShapeFunction(xi, a, val, dval_xi);

    // compute Jacobian:
    double phi[2];
    double J[2][2], JTinv[2][2];
    ComputeIsoparamMap(xi, nodal_coords, phi, J);

    // compute inverse of Jacobian transpose:
    InverseJacTranspose(J, JTinv);

    // compute determinant of Jacobian:
    detJ = Determinant(J);

    // compute gradient of shape function wrt x,y:
    dval[0] = JTinv[0][0] * dval_xi[0] + JTinv[0][1] * dval_xi[1];
    dval[1] = JTinv[1][0] * dval_xi[0] + JTinv[1][1] * dval_xi[1];
}   