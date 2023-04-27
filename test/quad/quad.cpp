

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <iomanip>

// (5 mins)
// Q1. Evaluate basis functions over the parametric element
void ComputeNhat(const double *xi, // in:  parametric coordinates over Khat
                 const int a,      // in:  basis function number = 0,1,2 or 3
                 double &Na,       // out: basis function value Na
                 double *dNa);     // out: gradient of the basis function grad_Na

// (10 mins)
// Q2. Compute the isoparametric map and its gradient
void Compute_isoparam_map(const double *xi,                        // in: parametric coordinates in Khat
                          const std::vector<double> &nodal_coords, // in: nodal coordinates in K
                          double phi[2],                           // out: isoparametric map
                          double grad_phi[2][2]);                  // out: gradient of the isoparametric map

// (15 mins)
// Q3. Compute the element force vector corresponding to the linear form F(v) = (f, v) with f(x,y) = x^2 + y^2
// use a 4 point quadrature rule
void ComputeFe(const std::vector<double> &nodal_coords, // in: nodal coordinates of the element x0,y0, x1,y1, x2,y2, x3,y3
               double Fe[4]);                           // out: element force vector

// (20 mins)
// Q4. Compute the element stiffness matrix corresponding to the bilinear form a(u,v) = (grad_u, grad_v)
// use a 4 point quadrature rule
void ComputeKe(const std::vector<double> &nodal_coords, // in: nodal coordinates of the element x0,y0, x1,y1, x2,y2, x3,y3
               double Ke[4][4]);                        // out: element stiffness matrix

int main()
{
  // test case 1. physical element = parametric element.
  std::vector<double> nodal_coords{-1., -1., -1., 1., 1., 1., -1., 1.};
  double Fe[4], Ke[4][4];
  ComputeFe(nodal_coords, Fe);
  ComputeKe(nodal_coords, Ke);
  std::cout << "Case 1: " << std::endl
            << "Fe: " << Fe[0] << " " << Fe[1] << " " << Fe[2] << " " << Fe[3] << std::endl;
  std::cout << "Ke: ";
  for (int i = 0; i < 4; ++i)
  {
    std::cout << std::endl;
    for (int j = 0; j < 4; ++j)
      std::cout << std::setprecision(4) << Ke[i][j] << " ";
  }
  std::cout << std::endl
            << std::endl;

  // test case 2
  nodal_coords = {0., 0., 1., 0., 2., 3., 1., 2.};
  ComputeFe(nodal_coords, Fe);
  ComputeKe(nodal_coords, Ke);
  std::cout << "Case 2: " << std::endl
            << "Fe: " << Fe[0] << " " << Fe[1] << " " << Fe[2] << " " << Fe[3] << std::endl;
  std::cout << "Ke: ";
  for (int i = 0; i < 4; ++i)
  {
    std::cout << std::endl;
    for (int j = 0; j < 4; ++j)
      std::cout << Ke[i][j] << " ";
  }
}

// Q1. Evaluate basis functions over the parametric element
void ComputeNhat(const double *xi, // in:  parametric coordinates over Khat
                 const int a,      // in:  basis function number = 0,1,2 or 3
                 double &Na,       // out: basis function value Na
                 double *dNa)      // out: gradient of the basis function grad_Na
{

  return;
}

// Q2. Compute the isoparametric map and its gradient
void Compute_isoparam_map(const double *xi,                        // in: parametric coordinates in Khat
                          const std::vector<double> &nodal_coords, // in: nodal coordinates in K
                          double phi[2],                           // out: isoparametric map
                          double grad_phi[2][2])                   // out: gradient of the isoparametric map
{




  return;
}

// Q3. Compute the element force vector
// use a 4 point quadrature rule
void ComputeFe(const std::vector<double> &nodal_coords, // in: nodal coordinates of the element x0,y0, x1,y1, x2,y2, x3,y3
               double Fe[4])                            // out: element force vector
{
  // ...
  return;
}

// Q4. Compute the element stiffness matrix
// use a 4 point quadrature rule
void ComputeKe(const std::vector<double> &nodal_coords, // in: nodal coordinates of the element x0,y0,x1,y1,...
               double Ke[4][4])                         // out: element stiffness matrix
{
  // ...
  return;
}
