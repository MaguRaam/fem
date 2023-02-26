// Solve Poisson equation with Dirchlet boundary condition on a triangular mesh:
#include "ex1.h"

// exact solution:
inline double uexact(const double &x, const double &y)
{
    return x + y;
}

// rhs function:
inline double rhs(const double &x, const double &y)
{
    return 0.0;
}

// compute area of triangle:
inline double triangle_area(const Vector3d &x, const Vector3d &y)
{

    double area = fabs(0.5 * (x(0) * (y(1) - y(2)) + x(1) * (y(2) - y(0)) + x(2) * (y(0) - y(1))));

    return area;
}

// compute element stiffness matrix:
inline Matrix3d compute_Klocal(const Vector3d &x, const Vector3d &y)
{

    double area = triangle_area(x, y);

    // derivative of nodal basis functions Nx and Ny:
    Vector3d Nx(y(1) - y(2), y(2) - y(0), y(0) - y(1));
    Vector3d Ny(x(2) - x(1), x(0) - x(2), x(1) - x(0));

    Nx = Nx / (2 * area);
    Ny = Ny / (2 * area);

    Matrix3d Klocal;
    Klocal = (Nx * Nx.transpose() + Ny * Ny.transpose()) * area;

    return Klocal;
}

// assemble stiffness matrix:
void assemble_K(SparseMatrix<double, RowMajor> &K, const MatrixX3i &conn, const VectorXd &x, const VectorXd &y)
{
    const int nElements = conn.rows();

    vector<Triplet<double>> tripleList;
    tripleList.reserve(nElements * 9);

    // loop over elements:
    for (int e = 0; e < nElements; ++e)
    {
        // get coordinates of the triangle nodes:
        Vector3d X, Y;
        for (int n = 0; n < 3; ++n)
        {
            X(n) = x(conn(e, n));
            Y(n) = y(conn(e, n));
        }

        // get local stiffness matrix given vertex coordinates of the element:
        Matrix3d Klocal = compute_Klocal(X, Y);

        // loop over nodes in the element:
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                tripleList.emplace_back(conn(e, i), conn(e, j), Klocal(i, j));
    }

    // build sparse matrix K:
    K.setFromTriplets(tripleList.begin(), tripleList.end());
}




int main()
{
    // load connectivity matrix:
    MatrixX3i conn;
    load_connectivity("mesh/ele.dat", conn);

    // load nodal coordinates:
    VectorXd x, y;
    load_coordinates("mesh/node.dat", x, y);

    // no of elements and nodes:
    const int nElements = conn.rows();
    const int nNodes = x.size();

    // compute exact solution:
    VectorXd Uexact(nNodes);
    for (int n = 0; n < nNodes; ++n)
        Uexact(n) = uexact(x(n), y(n));

    // assemble stiffness matrix K:
    SparseMatrix<double, RowMajor> K(nNodes, nNodes);
    assemble_K(K, conn, x, y);
    saveMarket(K, "K.mtx");

    // assemble force vector F:
    VectorXd F = VectorXd::Zero(nNodes);

    // set dirchlet boundary conditions:
    for (int n = 0; n < nNodes; ++n)
    {
        // check if node is located on the boundary:
        if (x(n) == 0.0 || x(n) == 1.0 || y(n) == 0.0 || y(n) == 1.0)
        {
            // make the nth row zeros:
            for (int j = 0; j < nNodes; ++j)
                K.coeffRef(n, j) = 0.0;

            // make the diagonal element 1:
            K.coeffRef(n, n) = 1.0;
        
            // set boundary conditon 
            F(n) = uexact(x(n), y(n));
        }
    }


    // solve KU = F:
    SparseLU<decltype(K)> solver;
    solver.analyzePattern(K);
    solver.factorize(K);
    VectorXd U = solver.solve(F);
    assert(solver.info() == Success);

    // compute abs Error:
    VectorXd Error = (U - Uexact).cwiseAbs();

    std::cout << "Linfty error = " << Error.maxCoeff() << std::endl;

    // write solution in .vtk:
    write_vtk(conn, x, y, U, Uexact, Error, "sol.vtk");

    return 0;
}
