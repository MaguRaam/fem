// pre and postprocessing tools
// compute L2 error

#pragma once

#include <Eigen/Dense>
#include "fe.h"

// compute L2 error:
template <typename UExact>
inline double L2error(UExact uexact, const Eigen::Vector3d &X, const Eigen::Vector3d &Y, const Eigen::Vector3d &U)
{
    // get quadrature pts and weights:
    const auto &nqpts = GaussQuadratureTriangle<3>::nqpts;
    const auto &wts = GaussQuadratureTriangle<3>::wts;
    const auto &r = GaussQuadratureTriangle<3>::r;
    const auto &s = GaussQuadratureTriangle<3>::s;

    double l2error = 0.0;

    // loop over quadrature pts:
    for (int q = 0; q < nqpts; ++q)
    {
        // get basis, basis gradients and Jacobian at quadpoint:
        auto [S, dSdx, dSdy, detJ] = isomap(r[q], s[q], X, Y);

        // map r, s -> x, y:
        double x = X.dot(S);
        double y = Y.dot(S);

        // finite element solution:
        double uh = U.dot(S);

        l2error += pow((uexact(x, y) - uh), 2.0 ) * detJ * 0.5 * wts[q];
    }

    return l2error;
}

//TODO test this function 
//?compute energy error:
template <typename GradUExact>
inline double energy_error(GradUExact grad_uexact, const Eigen::Vector3d &X, const Eigen::Vector3d &Y, const Eigen::Vector3d &U)
{
    // get quadrature pts and weights:
    const auto &nqpts = GaussQuadratureTriangle<3>::nqpts;
    const auto &wts = GaussQuadratureTriangle<3>::wts;
    const auto &r = GaussQuadratureTriangle<3>::r;
    const auto &s = GaussQuadratureTriangle<3>::s;

    double error = 0.0;

    // loop over quadrature pts:
    for (int q = 0; q < nqpts; ++q)
    {
        // get basis, basis gradients and Jacobian at quadpoint:
        auto [S, dSdx, dSdy, detJ] = isomap(r[q], s[q], X, Y);

        // map r, s -> x, y:
        double x = X.dot(S);
        double y = Y.dot(S);

        // finite element solution:
        Eigen::Vector2d grad_uh(U.dot(dSdx), U.dot(dSdy));

        error += (grad_uexact(x, y) - grad_uh).dot(grad_uexact(x, y) - grad_uh) * detJ * 0.5 * wts[q];
    }

    return error;
}








