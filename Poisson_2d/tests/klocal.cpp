// test local stiffness matrix:
#include "fe.h"
#include <iostream>

// compute element stiffness matrix:
inline Matrix3d compute_Klocal1(const Vector3d &X, const Vector3d &Y)
{
    Matrix3d Klocal = Matrix3d::Zero(); 

    // get quadrature pts and weights:
    const auto& nqpts = GaussQuadratureTriangle<1>::nqpts;
    const auto& wts = GaussQuadratureTriangle<1>::wts;
    const auto& r = GaussQuadratureTriangle<1>::r;
    const auto& s = GaussQuadratureTriangle<1>::s;

    // loop over quadrature pts and compute local stiffness matrix:
    for (int q = 0; q < nqpts; ++q)
    {
        // get gradients and Jacobian at quadpoint:
        auto [S, dSdx, dSdy, detJ] = isomap(r[q], s[q], X, Y);

        Klocal += (dSdx*dSdx.transpose() + dSdy*dSdy.transpose())*detJ*0.5*wts[q];  //! note we scale weight by 0.5:
    }

    return Klocal;
}


inline Matrix3d compute_Klocal2(const Vector3d &x, const Vector3d &y)
{
    // compute area of the triangle:
    double area = fabs(0.5 * (x(0) * (y(1) - y(2)) + x(1) * (y(2) - y(0)) + x(2) * (y(0) - y(1))));

    // derivative of nodal basis functions Nx and Ny:
    Vector3d Nx(y(1) - y(2), y(2) - y(0), y(0) - y(1));
    Vector3d Ny(x(2) - x(1), x(0) - x(2), x(1) - x(0));

    Nx = Nx / (2 * area);
    Ny = Ny / (2 * area);

	// dot(gradNi,gradNj)*area
    return (Nx * Nx.transpose() + Ny * Ny.transpose()) * area;
}


int main()
{
    // linear function:
    auto f = [](const double &x, const double &y)
    { return x + y; };

    Vector3d X{1.0, 2.0, 1.5}, Y{3.0, 3.5, 4.0};

    std::cout << (compute_Klocal1(X, Y) - compute_Klocal2(X, Y)) << std::endl;
}

/*
    map (r,s) -> (x, y)
    double x = X.dot(S);
    double y = Y.dot(S);
*/
