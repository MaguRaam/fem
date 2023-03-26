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
    const auto &nqpts = GaussQuadratureTriangle<1>::nqpts;
    const auto &wts = GaussQuadratureTriangle<1>::wts;
    const auto &r = GaussQuadratureTriangle<1>::r;
    const auto &s = GaussQuadratureTriangle<1>::s;

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

        l2error += pow((uexact(x, y) - uh), 2.0) * detJ * 0.5 * wts[q];
    }

    return l2error;
}

// TODO test this function
//?compute energy error:
template <typename GradUExact>
inline double energy_error(GradUExact grad_uexact, const Eigen::Vector3d &X, const Eigen::Vector3d &Y, const Eigen::Vector3d &U)
{
    // get quadrature pts and weights:
    const auto &nqpts = GaussQuadratureTriangle<1>::nqpts;
    const auto &wts = GaussQuadratureTriangle<1>::wts;
    const auto &r = GaussQuadratureTriangle<1>::r;
    const auto &s = GaussQuadratureTriangle<1>::s;

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

// Compute Energy functional of fem or interpolated solution:
template <typename F>
inline double energy_functional(const Eigen::Vector3d &X, const Eigen::Vector3d &Y, const Eigen::Vector3d &U, F f)
{
    // get quadrature pts and weights:
    const auto &nqpts = GaussQuadratureTriangle<1>::nqpts;
    const auto &wts = GaussQuadratureTriangle<1>::wts;
    const auto &r = GaussQuadratureTriangle<1>::r;
    const auto &s = GaussQuadratureTriangle<1>::s;

    double energy = 0.0;

    // loop over quadrature pts:
    for (int q = 0; q < nqpts; ++q)
    {
        // get basis, basis gradients and Jacobian at quadpoint:
        auto [S, dSdx, dSdy, detJ] = isomap(r[q], s[q], X, Y);

         // map r, s -> x, y:
        double x = X.dot(S);
        double y = Y.dot(S);

        // u and gradient of u at quadrature point:
        Eigen::Vector2d grad_u(U.dot(dSdx), U.dot(dSdy));
        double u = U.dot(S);

        energy += (0.5*grad_u.dot(grad_u) - f(x,y)*u)* detJ * 0.5 * wts[q];

    }

    return energy;
}

// Compute Energy functional of exact solution:
template <typename U, typename gradU, typename F>
inline double energy_functional(const Eigen::Vector3d &X, const Eigen::Vector3d &Y, U u, gradU grad_u, F f)
{
    // get quadrature pts and weights:
    const auto &nqpts = GaussQuadratureTriangle<1>::nqpts;
    const auto &wts = GaussQuadratureTriangle<1>::wts;
    const auto &r = GaussQuadratureTriangle<1>::r;
    const auto &s = GaussQuadratureTriangle<1>::s;

    double energy = 0.0;

    // loop over quadrature pts:
    for (int q = 0; q < nqpts; ++q)
    {
        // get basis, basis gradients and Jacobian at quadpoint:
        auto [S, dSdx, dSdy, detJ] = isomap(r[q], s[q], X, Y);

         // map r, s -> x, y:
        double x = X.dot(S);
        double y = Y.dot(S);

        energy += (0.5*grad_u(x, y).dot(grad_u(x, y)) - f(x,y)*u(x, y))* detJ * 0.5 * wts[q];

    }

    return energy;
}


