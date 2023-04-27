#ifndef FE_H
#define FE_H

#include <vector>
#include <cassert>

// four point quadratur formual for parametric square element:
static const int nqpts = 4;
static const double wts[nqpts] = {1.0, 1.0, 1.0, 1.0};
static const double xi[nqpts][2] = {{-0.5773502692, -0.5773502692},
                                    {0.5773502692, -0.5773502692},
                                    {0.5773502692, 0.5773502692},
                                    {-0.5773502692, 0.5773502692}};

// get shape function and derivatives over the parametric element:
void GetShapeFunction(const double *xi, // in:  parametric coordinates over Khat
                      const int a,      // in:  basis function number = 0,1,2 or 3
                      double &Na,       // out: basis function value Na
                      double *dNa);     // out: gradient of the basis function grad_Na

// compute iso parametric map and Jacobian of the map:
void ComputeIsoparamMap(const double *xi,                        // in: parametric coordinates in Khat
                        const std::vector<double> &nodal_coords, // in: nodal coordinates in K
                        double phi[2],                           // out: isoparametric map
                        double grad_phi[2][2]);                  // out: gradient of the isoparametric map

// compute determinant of Jacobian matrix:
double Determinant(const double grad_phi[2][2] // gradient of isoparametric map
);

// inverse of Jacobian transpose:
void InverseJacTranspose(const double J[2][2], // in:  Jacobian of map
                         double JTinv[2][2]    // out: inverse of Jacobian transpose

);

// shape functions and derivatives over the physical element:
void GetShapeFunction(const double *xi,                        // in:  parametric coordinates over Khat
                      int a,                                   // in:  basis function number = 0,1,2 or 3
                      const std::vector<double> &nodal_coords, // in: nodal coordinates in K
                      double &val,                             // out: basis function value
                      double *dval,                            // out: gradient of the basis function over a physical element
                      double &detJ                             // determinant of Jacobian
);

#endif