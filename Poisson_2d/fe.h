// Quadrature formula for triangle element:
// Linear basis functions on triangle element P1(K):
// Isoparametric map:
// Integrate arbirtary function over an element K:

#pragma once

#include <array>
#include <Eigen/Dense>

// Gauss points and weights on a reference triangle:
// give polynomial degree as template parameter:
// !Note: note the weights are scaled so they sum to one, so divide det(J) by 2:

template <int degree>
struct GaussQuadratureTriangle;

template <>
struct GaussQuadratureTriangle<1>
{
    static constexpr int nqpts = 1;
    static constexpr std::array<double, nqpts> wts{1.0};
    static constexpr std::array<double, nqpts> r{1.0 / 3.0};
    static constexpr std::array<double, nqpts> s{1.0 / 3.0};
};

template <>
struct GaussQuadratureTriangle<2>
{
    static constexpr int nqpts = 3;
    static constexpr std::array<double, nqpts> wts{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    static constexpr std::array<double, nqpts> r{1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0};
    static constexpr std::array<double, nqpts> s{1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0};
};

template <>
struct GaussQuadratureTriangle<3>
{
    static constexpr int nqpts = 4;
    static constexpr std::array<double, nqpts> wts{-27.0 / 48.0, 25.0 / 48.0, 25.0 / 48.0, 25.0 / 48.0};
    static constexpr std::array<double, nqpts> r{1.0 / 3.0, 0.2, 0.6, 0.2};
    static constexpr std::array<double, nqpts> s{1.0 / 3.0, 0.2, 0.2, 0.6};
};

// given r,s in reference triangle element, returns basis function and gradient:
//(r, s) -> {S, dSdr, dSds}
inline auto basis_function(const double &r, const double &s) -> std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>
{
    return {{1 - r - s, r, s}, {-1.0, 1.0, 0.0}, {-1.0, 0.0, 1.0}};
}

// given r,s in reference triangle element and nodal coordinates X,Y of physical triangle, returns basis function and gradient wrt x,y and determinant of Jacogian:
// (r, s) -> {S, dSdx, dSdy, detJ}
inline auto isomap(const double &r, const double &s, const Eigen::Vector3d &X, const Eigen::Vector3d &Y) -> std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, double>
{
    auto [S, dSdr, dSds] = basis_function(r, s);

    // compute Jacobian matrix:
    Eigen::Matrix2d J;
    J << X.dot(dSdr), Y.dot(dSdr),
        X.dot(dSds), Y.dot(dSds);

    // compute determinant of Jacobian matrix:
    double detJ = J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);

    // compute gradient of basis functions dSdx and dSdy:
    Eigen::Vector3d dSdx = (J(1, 1) * dSdr - J(0, 1) * dSds) / detJ;
    Eigen::Vector3d dSdy = (-J(1, 0) * dSdr + J(0, 0) * dSds) / detJ;

    return {S, dSdx, dSdy, detJ};
}

