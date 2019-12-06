#include <vector>
#include <timestepping.h>
#include <assemble.h>
#include <eigen3/Eigen/Dense>
enum save_vtk
{
    YES,
    NO
};

void StepRK4(Eigen::VectorXd &u_init, double curr_time, double dt, int N_DOF, Assemble& assemble)
{
    Eigen::VectorXd k1(N_DOF), k2(N_DOF), k3(N_DOF), k4(N_DOF);

    assemble.GradFunctionForRK(u_init, curr_time, k1);
    assemble.GradFunctionForRK(u_init + 0.5 * dt * k1, curr_time + 0.5 * dt, k2);
    assemble.GradFunctionForRK(u_init + 0.5 * dt * k2, curr_time + 0.5 * dt, k3);
    assemble.GradFunctionForRK(u_init + dt * k3, curr_time + dt, k4);

    u_init += u_init + (1.0 / 6.0) * dt * (k1 + 2 * k2 + 2 * k3 + k4);
}

/* 
    Eigen::VectorXd du_dt(N_DOF), k1(N_DOF), k2(N_DOF), k3(N_DOF), k4(N_DOF);
    // rhs.resize(N_DOF);
    // u.resize(N_DOF);
    // k1.resize(N_DOF);
    // k2.resize(N_DOF);
    // k3.resize(N_DOF);
    // k4.resize(N_DOF);
    // du_dt.resize(N_DOF);

    Eigen::VectorXd u_k1(N_DOF), u_k2(N_DOF), u_k3(N_DOF);

    AssembleRHS(rhs, u);
    ApplyInlet(rhs, u, val[0]);
    Getdu_dtForRK(u, rhs, k1);
    u_k1 = u + (0.5 * dt * k1);

    AssembleRHS(rhs, u_k1);
    ApplyInlet(rhs, u_k1, val[1]);
    Getdu_dtForRK(u_k1, rhs, k2);
    u_k2 = u + (0.5 * dt * k2);

    AssembleRHS(rhs, u_k2);
    ApplyInlet(rhs, u_k2, val[1]);
    Getdu_dtForRK(u_k2, rhs, k3);
    u_k3 = u + (dt * k3);

    AssembleRHS(rhs, u_k3);
    ApplyInlet(rhs, u_k3, val[2]);
    Getdu_dtForRK(u_k2, rhs, k4);
    u = (u + (1.0 / 6.0) * dt * (k1 + 2 * k2 + 2 * k3 + k4));
*/