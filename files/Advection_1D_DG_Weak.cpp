#include <iostream>
#include <basis_function.h>
#include <list>
#include <fstream>
#include <matrix.h>
#include <mesh.h>
#include <assemble.h>
#include <assembleDG.h>
#include <iomanip> 
#include <solver.h>
#include <post_process.h>
#include <cmath>
#include <limits>
#include <algorithm>
#include <eigen3/Eigen/Dense>
#define DEBUG_PRINT_EIGEN_VEC(u)                     \
    for (int i = 0; i < order + 1; i++)              \
    {                                                \
        for (int j = 0; j < mesh.get_N_Cells(); j++) \
        {                                            \
            cout << u[i + j * N_Local_DOF] << "\t";  \
        }                                            \
        cout << endl;                                \
    }
using namespace std;

int main(int argc, char **argv){
    int order = stoi(argv[1]);
    int N_Local_DOF = order + 1;
    Mesh mesh(0, 2, stoi(argv[2]), order);
    int N_DOF = N_Local_DOF * mesh.get_N_Cells();
    
    Eigen::VectorXd u, rhs, du_dt, k1, k2, k3, k4; 
    LegendreBasis Lbasis(order);
    vector<double>& Lp = Lbasis.LGLPoints;
    AssembleDG assemble(Lbasis, order, mesh, 0.0);

    assemble.AssembleMkS();

    rhs.resize(N_DOF);
    u.resize(N_DOF);
    k1.resize(N_DOF);
    k2.resize(N_DOF);
    k3.resize(N_DOF);
    k4.resize(N_DOF);
    du_dt.resize(N_DOF);

    fill(rhs.data(), rhs.data() + N_DOF, 0.0);
    fill(u.data(), u.data() + N_DOF, 0.0);

    double x1 = mesh.get_cell(0).get_x1(), x2 = mesh.get_cell(0).get_x2();
    auto mapping = [x1, x2](double si) {
        return ((1 - si) / 2) * x1 + ((1 + si) / 2) * x2;
    };

    double cfl = 0.75,
    dx = numeric_limits<double>::max();
    for (int i = 0; i < order; i++)
        dx = min(dx, mapping(Lbasis.LGLPoints[i + 1]) - mapping(Lbasis.LGLPoints[i]));

    double t = 0.0, dt = (cfl / (2 * M_PI)) * dx / 2; // Careful with the step size
    int timesteps = stoi(argv[3]);

    vector<double> X;
    for (int i = 0; i < mesh.get_N_Cells(); i++){
        X.push_back(mesh.get_cell(i).get_x1());
        X.push_back(mesh.get_cell(i).get_x2());
    }
    X.reserve(mesh.get_N_Cells() * N_Local_DOF);
    

    for (int t_i = 0; t_i < timesteps; t_i++)
    {
        double val[3];
        val[0] = -1.0 * sin(2 * M_PI * t);
        val[1] = -1.0 * sin(2 * M_PI * (t + 0.5 * dt));
        val[2] = -1.0 * sin(2 * M_PI * (t + dt));

        assemble.StepRKWeak(u, rhs, t, dt, val);

        vector<double> Y;
        Y.reserve(mesh.get_N_Cells() * N_Local_DOF);
        for (int i = 0; i < mesh.get_N_Cells(); i++)
        {
            Y.push_back(u(i * N_Local_DOF));
            Y.push_back(u((i + 1) * N_Local_DOF - 1));
        }
        stringstream ss;
        ss << "test-" << t_i + 1;
        write_vtk(X, Y, ss.str().c_str());
        ss.clear();
        t += dt;
    }
    cout << setprecision(3);
    for (int i = 0; i < order + 1; i++)
    {
        for (int j = 0; j < mesh.get_N_Cells(); j++){
            cout << u[i + j * N_Local_DOF] << "\t";
        }
        cout << endl;
    }

    

    double each_cell_disc = 10;
    basis_function lagrange(order, Lp);
    
    return 0;
}
