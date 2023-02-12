#include "L2.h"

// Gauss points and weights on a reference triangle:
//! Compile time computation of quadrature pts and weights:
template <int npts>
struct GaussQuadrature;

// 1 pt Gauss-Quadrature on a triangle:
template <>
struct GaussQuadrature<1>
{
    static constexpr std::array<double, 1> wts(){
        return {0.5};
    }

    static constexpr std::array<double, 1> r(){
        return {1.0 / 3.0};
    }

    static constexpr std::array<double, 1> s(){
        return {1.0 / 3.0};
    }

};

// 3 pt Gauss-Quadrature on a triangle:
template <>
struct GaussQuadrature<3>
{
    static constexpr std::array<double, 3> wts(){
        return {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0};
    }

    static constexpr std::array<double, 3> r(){
        return {1.0/6.0, 2.0/3.0, 1.0/6.0};
    }

    static constexpr std::array<double, 3> s(){
        return {1.0/6.0, 1.0/6.0, 2.0/3.0};
    }

};

// returns shape functions and its derivatives wrt r,s on a reference triangle:
inline auto shape_function(double r, double s) -> std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>
{
    Eigen::Vector3d S{1 - r - s, r, s}, dSdr{-1.0, 1.0, 0.0}, dSds{-1.0, 0.0, 1.0};
    return {S, dSdr, dSds};
}

// compute shape functions, derivatives wrt x, y and detJ at r,s given vertices of physical triangle:
inline auto isoparametric(double r, double s, const Eigen::Vector3d &x, const Eigen::Vector3d &y) -> std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, double>
{
    // get shape function and its derivatives wrt r,s given r,s on a reference triangle:
    auto [S, dSdr, dSds] = shape_function(r, s);

    // Compute Jacobian matrix:
    Eigen::Matrix2d J;

    J(0, 0) = dSdr.transpose() * x;
    J(0, 1) = dSdr.transpose() * y;
    J(1, 0) = dSds.transpose() * x;
    J(1, 1) = dSds.transpose() * y;

    // Compute determinant of Jacobian:
    double detJ = J.determinant();

    // Compute derivative of shape functions wrt x, y at r, s:
    Eigen::Vector3d dSdx, dSdy;
    dSdx = (J(1, 1) * dSdr - J(0, 1) * dSds) / detJ;
    dSdy = (-J(1, 0) * dSdr + J(0, 0) * dSds) / detJ;

    return {S, dSdx, dSdy, detJ};
}

// compute local mass matrix:
inline Eigen::Matrix3d mass_matrix(const Eigen::Vector3d &x, const Eigen::Vector3d &y)
{
    Eigen::Matrix3d Mlocal = Eigen::Matrix3d::Zero();
    
    // get quad points and weights on a reference triangle:
    constexpr int nqpts = 3;
    auto wt = GaussQuadrature<nqpts>::wts();
    auto r = GaussQuadrature<nqpts>::r();
    auto s = GaussQuadrature<nqpts>::s();

    // loop over quadrature pts:
    for (int q = 0; q < nqpts; ++q){

        auto [S, dSdx, dSdy, detJ] = isoparametric(r[q], s[q], x, y);

        Mlocal += S*S.transpose()*detJ*wt[q];

    }
    return Mlocal;
}




int main()
{
    std::vector<point> coordinates;
    std::vector<std::array<int, 3>> connectivity;

    // read mesh:
    read_mesh("mesh/coordinates-0.dat", coordinates, "mesh/connectivity-0.dat", connectivity);

    // write mesh:
    write_mesh("mesh.dat", coordinates, connectivity);

    // coordinates of a triangle:
    Eigen::Vector3d x{0.0, 1.0, 0.0}, y{0.0, 0.0, 1.0};

    std::cout << mass_matrix(x, y) << std::endl;


    return 0;
}
