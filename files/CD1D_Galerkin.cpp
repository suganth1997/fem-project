#include <iostream>
#include <basis_function.h>
#include <list>
#include <fstream>
#include <matrix.h>
#include <mesh.h>
#include <assemble.h>
#include <iomanip> 
#include <solver.h>
#include <post_process.h>
#include <cmath>
#include <limits>
#include <algorithm>
#include <functional>
#include <timestepping.h>
#include <eigen3/Eigen/Dense>

using namespace std;


int main(int argc, char **argv)
{
    int N_Cells = 10, order = 2, N_Local_DOF = order + 1;
    double l = 1.0;
    Mesh mesh(0, l, N_Cells, order);
    int N_DOF = N_DOF = mesh.get_N_DOF();
    sparsity_pattern pattern(mesh, N_Local_DOF);
    pattern.compute_sparsity();

    Matrix matrix(mesh.get_N_DOF(), N_Local_DOF, pattern), mass_matrix(mesh.get_N_DOF(), N_Local_DOF, pattern);

    vector<double> X;
    X.reserve(mesh.get_N_DOF());

    Assemble assemble(mesh, matrix, mass_matrix, order);
    assemble.AssembleGlobalSystemCD1D(stof(argv[1]), stof(argv[2]));
    assemble.AssembleMassMatrix();

    Eigen::VectorXd sol(mesh.get_N_DOF());
    
    sol.resize(mesh.get_N_DOF());

    double val = 100; /* Neumann boundary at the right */

    for (int i = 0; i < mesh.get_N_Cells(); i++){
        const cell_1d &cell = mesh.get_cell(i);
        const vector<double> &vertices = cell.get_vertices();
        int begin = mesh.global_DOF(i);
        for (int j = 0; j < order + 1; j++)
        {
            sol(begin + j) = vertices[j] * val / l;
        }
    }

    Eigen::VectorXd rhs;
    rhs.resize(N_DOF);

    const vector<int> &RowPtr = mass_matrix.get_RowPtr(), &colIndex = mass_matrix.get_colIndex();

    Eigen::MatrixXd K(N_DOF - 2, N_DOF - 2), M(N_DOF - 2, N_DOF - 2), M_inv;
    cout << setprecision(3);
    for (int i = 1; i < N_DOF - 1; i++){
        for (int j = RowPtr[i]; j < RowPtr[i + 1]; j++){
            if(colIndex[j] == 0 or colIndex[j] == N_DOF - 1)
                continue;
            K(i - 1, colIndex[j] - 1) = matrix.entries[j];
            M(i - 1, colIndex[j] - 1) = mass_matrix.entries[j];
        }
    }

    double dt = 0.025;
    int n_steps = 200;

    Eigen::VectorXd rhs_diff(N_DOF - 2);
    
    sol(N_DOF - 1) = val;

    M_inv = M.inverse();

    for (int i = 0; i < mesh.get_N_DOF(); i++)
    {
        X.push_back(i * l / (mesh.get_N_DOF() - 1));
    }

    // write_vtk_line(X, sol, "test-0");

    vector<double> x;
    for (int i = 0; i < mesh.get_N_Cells(); i++){
        x.push_back(mesh.get_cell(i).get_x1());
        x.push_back(mesh.get_cell(i).get_x2());
    }

    for (int i = 0; i < n_steps; i++)
    {
        stringstream ss;
        ss << "test-" << i + 1;

        for (int i = 0; i < N_DOF - 2; i++){
            rhs_diff[i] = (mass_matrix.get(i + 1, N_DOF - 1) + dt * matrix.get(i + 1, N_DOF - 1)) * val;
        }

        Eigen::VectorXd sol_unknown(N_DOF - 2);
        memcpy(sol_unknown.data(), sol.data() + 1, (N_DOF - 2) * sizeof(double));

        sol_unknown = (M + (dt * K)).inverse() * (M * sol_unknown - rhs_diff);

        memcpy(sol.data() + 1, sol_unknown.data(), (N_DOF - 2) * sizeof(double));

        // write_vtk_line(X, sol, ss.str().c_str());
    }

    gnuplot(X, sol);
    return 0;
}