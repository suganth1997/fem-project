#include <assembleDG.h>
#include <assemble.h>
#include <vector>

using namespace std;

AssembleDG::AssembleDG(LegendreBasis &Lbasis, int order, const Mesh& mesh, double alpha) : 
Lbasis(Lbasis), order(order), N_Local_DOF(order + 1), mesh(mesh), alpha(alpha), N_DOF((order+1)*mesh.get_N_Cells()) {

}

void AssembleDG::AssembleMkSLagrange(basis_function &lagrange){
    const vector<double> &xi = lagrange.get_nodal_X();

    M = Eigen::MatrixXd::Constant(N_Local_DOF, N_Local_DOF, 0.0);
    S = Eigen::MatrixXd::Constant(N_Local_DOF, N_Local_DOF, 0.0);

    QuadPointData1D qpdata(order + 1, order);

    vector<pair<double, pair<vector<double>, vector<double>>>>& quadWtsFnsGrads = qpdata.quadratureWeightsFunctionsAndGradients;

    for (int w_i = 0; w_i < qpdata.n_quad_pts; w_i++){
        double quad_wt = quadWtsFnsGrads[w_i].first;
        vector<double> &shape_fns = quadWtsFnsGrads[w_i].second.first;
        vector<double> &shape_grads = quadWtsFnsGrads[w_i].second.second;

        for (int i = 0; i < N_Local_DOF; i++)
        {
            for (int j = 0; j < N_Local_DOF; j++){
                M(i, j) += quad_wt * shape_fns[i] * shape_fns[j];
                S(i, j) += quad_wt * shape_fns[i] * shape_grads[j]; //j == 0 or j == N_Local_DOF - 1 ? fabs(j - (N_Local_DOF - 1)) : j
            }
        }
    }

    Dr = M.inverse() * S;
    // for (int i = 0; i < N_Local_DOF; i++){
    //     for (int j = 0; j < N_Local_DOF; j++){
    //         cout << S(i, j) << "\t";
    //     }
    //     cout << endl;
    // }
    // exit(0);
    // S = -1.0 * S;
}

void AssembleDG::AssembleMkS(){
    vector<double> &r = Lbasis.LGLPoints;

    Eigen::MatrixXd V(N_Local_DOF, N_Local_DOF), Vr(N_Local_DOF, N_Local_DOF);

    for (int i = 0; i < N_Local_DOF; i++)
    {
        for (int j = 0; j < N_Local_DOF; j++){
            V(i, j) = Lbasis.Pn(r[i], j, 0, 0);
            Vr(i, j) = Lbasis.GradPn(r[i], j, 0, 0);
        }
    }

    M = (V * V.transpose()).inverse();
    Dr = Vr * V.inverse();
    S = M * Dr;

    hk.resize(mesh.get_N_Cells(), 0.0);

    for (int i = 0; i < mesh.get_N_Cells(); i++){
        const cell_1d &cell = mesh.get_cell(i);
        hk[i] = (cell.get_x2() - cell.get_x1()) / 2.0;
    }
}

void AssembleDG::AssembleRHS(Eigen::VectorXd& rhs, const Eigen::VectorXd& u){
    int N_Cells = mesh.get_N_Cells();
    double u_star;
    for (int i = 0; i < N_Cells; i++)
    {
        if(!i){
            u_star = ((u[N_Local_DOF - 1] + u[N_Local_DOF]) / 2) + ((1 - alpha) * (u[N_Local_DOF - 1] - u[N_Local_DOF]) / 2);
            rhs[N_Local_DOF - 1] = 2 * M_PI * (u[N_Local_DOF - 1] - u_star);
        }
        else if(i == N_Cells - 1){
            u_star = ((u[(N_Local_DOF * i) - 1] + u[N_Local_DOF * i]) / 2) +
                     ((1 - alpha) * (u[(N_Local_DOF * i) - 1] - u[N_Local_DOF * i]) / 2);

            rhs[N_Local_DOF * i] = -2 * M_PI * (u[N_Local_DOF * i] - u_star);
        }
        else{
            u_star = ((u[(N_Local_DOF * i) - 1] + u[N_Local_DOF * i]) / 2) +
                     ((1 - alpha) * (u[(N_Local_DOF * i) - 1] - u[N_Local_DOF * i]) / 2);

            rhs[N_Local_DOF * i] = -2 * M_PI * (u[N_Local_DOF * i] - u_star);

            u_star = ((u[N_Local_DOF * (i + 1) - 1] + u[N_Local_DOF * (i + 1)]) / 2) +
                     ((1 - alpha) * (u[N_Local_DOF * (i + 1) - 1] - u[N_Local_DOF * (i + 1)]) / 2);

            rhs[N_Local_DOF * (i + 1) - 1] = 2 * M_PI * (u[N_Local_DOF * (i + 1) - 1] - u_star);
        }
    }
    // rhs[0] = rhs[N_Local_DOF * N_Cells - 1] = 0.0;
}

void AssembleDG::ApplyInlet(Eigen::VectorXd& rhs, const Eigen::VectorXd& u, double val){
    double u_star = (val + u[0]) / 2 + ((1 - alpha) * (val - u[0]) / 2);

    rhs[0] = -2 * M_PI * (u[0] - val); // val instead of u_star stopped divergence, dont know why!!
    rhs[N_Local_DOF * mesh.get_N_Cells() - 1] = 0.0;
}

void AssembleDG::Getdu_dtForRK(const Eigen::VectorXd& u, const Eigen::VectorXd& rhs, Eigen::VectorXd &du_dt){
    Eigen::MatrixXd M_inv = M.inverse();
    Eigen::VectorXd uk(order + 1), rk(order + 1);

    int N_Cells = mesh.get_N_Cells();
    for (int i = 0; i < N_Cells; i++){
        const cell_1d &cell = mesh.get_cell(i);
        double hk = cell.get_x2() - cell.get_x1();

        // memcpy(uk.data(), u.data() + (i * N_Local_DOF), N_Local_DOF * sizeof(double));
        for (int j = 0; j < N_Local_DOF; j++)
        {
            uk[j] = u[i * N_Local_DOF + j];
            rk[j] = rhs[i * N_Local_DOF + j];
        }

        Eigen::VectorXd dudt = (2 / hk) * ((M_inv * rk) - (2 * M_PI * (Dr * uk)));

        for (int j = 0; j < N_Local_DOF; j++)
            du_dt[i * N_Local_DOF + j] = dudt[j];
    }
}

void AssembleDG::StepRK(Eigen::VectorXd &u, Eigen::VectorXd &rhs, double t, double dt, double val[]){
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
    
}


void AssembleDG::AssembleRHSWeak(Eigen::VectorXd& rhs, const Eigen::VectorXd& u, double coeff, double rhs_val){
    int N_Cells = mesh.get_N_Cells();
    double u_star;
    for (int i = 0; i < N_Cells; i++)
    {
        if(!i){
            u_star = ((u[N_Local_DOF - 1] + u[N_Local_DOF]) / 2) + ((1 - alpha) * (u[N_Local_DOF - 1] - u[N_Local_DOF]) / 2);
            rhs[N_Local_DOF - 1] = -coeff * u_star;
        }
        else if(i == N_Cells - 1){
            u_star = ((u[(N_Local_DOF * i) - 1] + u[N_Local_DOF * i]) / 2) +
                     ((1 - alpha) * (u[(N_Local_DOF * i) - 1] - u[N_Local_DOF * i]) / 2);

            rhs[N_Local_DOF * i] = coeff * u_star;
        }
        else{
            u_star = ((u[(N_Local_DOF * i) - 1] + u[N_Local_DOF * i]) / 2) +
                     ((1 - alpha) * (u[(N_Local_DOF * i) - 1] - u[N_Local_DOF * i]) / 2);

            rhs[N_Local_DOF * i] = coeff * u_star;

            u_star = ((u[N_Local_DOF * (i + 1) - 1] + u[N_Local_DOF * (i + 1)]) / 2) +
                     ((1 - alpha) * (u[N_Local_DOF * (i + 1) - 1] - u[N_Local_DOF * (i + 1)]) / 2);

            rhs[N_Local_DOF * (i + 1) - 1] = -coeff * u_star;
        }
    }
}

void AssembleDG::ApplyInletWeak(Eigen::VectorXd& rhs, const Eigen::VectorXd& u, double val, double coeff){
    // double u_star = ((val + u[0]) / 2) + ((1 - alpha) * (val - u[0]) / 2);

    // rhs[0] = 0;
    // //coeff *u_star; // val instead of u_star stopped divergence, dont know why!!
    // rhs[N_Local_DOF * mesh.get_N_Cells() - 1] = -coeff * u[u.size() - 1];
    double u_star = ((val + u[0]) / 2) + ((1 - alpha) * (val - u[0]) / 2);

    rhs[0] = 2.0 * M_PI * u_star; // val instead of u_star stopped divergence, dont know why!!
    rhs[N_Local_DOF * mesh.get_N_Cells() - 1] = -2.0 * M_PI * u[u.size() - 1];
}

void AssembleDG::ApplyOutletWeak(Eigen::VectorXd& rhs, const Eigen::VectorXd& u, double val, double coeff){
    double u_star = ((val + u[u.size() - 1]) / 2) + ((1 - alpha) * (u[u.size() - 1] - val) / 2);

    // rhs[0] = 2.0 * M_PI * u_star; // val instead of u_star stopped divergence, dont know why!!
    rhs[N_Local_DOF * mesh.get_N_Cells() - 1] = 0;
    //-coeff *u_star;
}

void AssembleDG::Getdu_dtForRKWeak(const Eigen::VectorXd &u, const Eigen::VectorXd &rhs, Eigen::VectorXd &du_dt){
    // cout << M << endl
    //      << endl
    //      << S << endl;
    // exit(0);
    Eigen::MatrixXd M_inv = M.inverse();
    Eigen::MatrixXd M_invS = (M_inv * S.transpose());
    Eigen::VectorXd uk(order + 1), rk(order + 1);

    int N_Cells = mesh.get_N_Cells();
    for (int i = 0; i < N_Cells; i++){
        const cell_1d &cell = mesh.get_cell(i);
        double hk = cell.get_x2() - cell.get_x1();

        // memcpy(uk.data(), u.data() + (i * N_Local_DOF), N_Local_DOF * sizeof(double));
        for (int j = 0; j < N_Local_DOF; j++)
        {
            uk[j] = u[i * N_Local_DOF + j];
            rk[j] = rhs[i * N_Local_DOF + j];
        }

        Eigen::VectorXd dudt = (2 / hk) * ((M_inv * rk) + (2 * M_PI * (M_invS * uk)));

        for (int j = 0; j < N_Local_DOF; j++)
            du_dt[i * N_Local_DOF + j] = dudt[j];
    }
}

void AssembleDG::StepRKWeak(Eigen::VectorXd &u, Eigen::VectorXd &rhs, double t, double dt, double val[]){
    Eigen::VectorXd du_dt(N_DOF), k1(N_DOF), k2(N_DOF), k3(N_DOF), k4(N_DOF);
    // rhs.resize(N_DOF);
    // u.resize(N_DOF);
    // k1.resize(N_DOF);
    // k2.resize(N_DOF);
    // k3.resize(N_DOF);
    // k4.resize(N_DOF);
    // du_dt.resize(N_DOF);

    Eigen::VectorXd u_k1(N_DOF), u_k2(N_DOF), u_k3(N_DOF);

    AssembleRHSWeak(rhs, u, 2 * M_PI);
    ApplyInletWeak(rhs, u, val[0], 2 * M_PI);
    Getdu_dtForRKWeak(u, rhs, k1);
    u_k1 = u + (0.5 * dt * k1);

    AssembleRHSWeak(rhs, u_k1, 2 * M_PI);
    ApplyInletWeak(rhs, u_k1, val[1], 2 * M_PI);
    Getdu_dtForRKWeak(u_k1, rhs, k2);
    u_k2 = u + (0.5 * dt * k2);

    AssembleRHSWeak(rhs, u_k2, 2 * M_PI);
    ApplyInletWeak(rhs, u_k2, val[1], 2 * M_PI);
    Getdu_dtForRKWeak(u_k2, rhs, k3);
    u_k3 = u + (dt * k3);

    AssembleRHSWeak(rhs, u_k3, 2 * M_PI);
    ApplyInletWeak(rhs, u_k3, val[2], 2 * M_PI);
    Getdu_dtForRKWeak(u_k2, rhs, k4);
    u = (u + (1.0 / 6.0) * dt * (k1 + 2 * k2 + 2 * k3 + k4));
}