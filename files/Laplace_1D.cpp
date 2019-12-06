#include <iostream>
#include <basis_function.h>
#include <list>
#include <fstream>
#include <matrix.h>
#include <mesh.h>
#include <assemble.h>
#include <solver.h>
#include <post_process.h>
#include <eigen3/Eigen/Sparse>
using namespace std;

int main(int argc, char **argv)
{
    int N_Cells = 20, order = 2, N_Local_DOF = order + 1;
    double l = 1.0;
    Mesh mesh(0.0, l, N_Cells, order);

    sparsity_pattern pattern(mesh, N_Local_DOF);
    pattern.compute_sparsity();

    Matrix matrix(mesh.get_N_DOF(), N_Local_DOF, pattern);

    Assemble assemble(mesh, matrix, matrix, order);
    assemble.AssembleGlobalSystem(1.0);

    vector<double> rhs, sol;
    rhs.resize(mesh.get_N_DOF(), 0.0);
    sol.resize(mesh.get_N_DOF(), 0.0);
    assemble.ApplyBoundaryConditionSymmetric(LEFT, DIRICHLET, 0.0, rhs);
    assemble.ApplyBoundaryConditionSymmetric(RIGHT, DIRICHLET, 1.0, rhs);

    sol[mesh.get_N_DOF() - 1] = 1.0;
    matrix.print();

    Solver solver(10000, 1e-10);
    solver.solve_cg(matrix, rhs, sol, mesh.get_N_DOF());

    vector<double> X;
    X.reserve(mesh.get_N_DOF());
    for (int i = 0; i < mesh.get_N_DOF(); i++)
    {
        X.push_back(i * l / (mesh.get_N_DOF() - 1));
    }

    gnuplot(X, sol);

    // write_vtk_line(X, sol, "laplace-1d");
    return 0;
}