#ifndef FE_H
#define FE_H

// quadarture pts and weights:
static const int nqpts = 1;
static const double wts[nqpts] = {1.0};
static const double xi[nqpts][2] = {{1.0 / 3.0, 1.0 / 3.0}};

// shape functions and derivatives over the parametric element:
void GetShapeFunction(const double *xi, int a, double &val, double *dval);

// shape functions and derivatives over the physical element:
void GetShapeFunction(const double *xi, int a, const double *nodal_coords, double &val, double *dval, double &detJ);

// map parametric element to physical element:
void isomap(const double *xi, double *x, const double *nodal_coords);


#endif