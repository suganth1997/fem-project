#include <iostream>
#include <assemble.h>
#include <basis_function.h>
#include <stdexcept>
#include <cmath> 
#include <assert.h>
#include <solver.h>
#include <eigen3/Eigen/Dense>
using namespace std;

void QuadPointData1D::GetWeightsAndPoints()
{
    switch(n_quad_pts){
        case 1:{
            quadWeightsAndPoints.reserve(1);
            quadWeightsAndPoints.emplace_back(2.0, 0.0);
            break;
        }
        
        case 2:{
            quadWeightsAndPoints.reserve(2);
            quadWeightsAndPoints.emplace_back(1.0, -0.5773502691896257);
            quadWeightsAndPoints.emplace_back(1.0, 0.5773502691896257);
            break;
        }

        case 3:{
            quadWeightsAndPoints.reserve(3);
            quadWeightsAndPoints.emplace_back(0.5555555555555556, -0.7745966692414834);
            quadWeightsAndPoints.emplace_back(0.8888888888888888, 0.0);
            quadWeightsAndPoints.emplace_back(0.5555555555555556, 0.7745966692414834);
            break;
        }

        case 4:{
            quadWeightsAndPoints.reserve(4);
            quadWeightsAndPoints.emplace_back(0.3478548451374538, -0.8611363115940526);
            quadWeightsAndPoints.emplace_back(0.6521451548625461, -0.3399810435848563);
            quadWeightsAndPoints.emplace_back(0.6521451548625461, 0.3399810435848563);
            quadWeightsAndPoints.emplace_back(0.3478548451374538, 0.8611363115940526);
            break;
        }

        case 5:{
            quadWeightsAndPoints.reserve(5);
            quadWeightsAndPoints.emplace_back(0.2369268850561891, -0.906179845938664);
            quadWeightsAndPoints.emplace_back(0.4786286704993665, -0.5384693101056831);
            quadWeightsAndPoints.emplace_back(0.5688888888888889, 0.0);
            quadWeightsAndPoints.emplace_back(0.4786286704993665, 0.5384693101056831);
            quadWeightsAndPoints.emplace_back(0.2369268850561891, 0.906179845938664);
            break;
        }

        case 6:{
            quadWeightsAndPoints.reserve(6);
            quadWeightsAndPoints.emplace_back(0.1713244923791704, -0.932469514203152);
            quadWeightsAndPoints.emplace_back(0.3607615730481386, -0.6612093864662645);
            quadWeightsAndPoints.emplace_back(0.467913934572691, -0.2386191860831969);
            quadWeightsAndPoints.emplace_back(0.467913934572691, 0.2386191860831969);
            quadWeightsAndPoints.emplace_back(0.3607615730481386, 0.6612093864662645);
            quadWeightsAndPoints.emplace_back(0.1713244923791704, 0.932469514203152);
        }

        case 7:{
            quadWeightsAndPoints.reserve(7);
            quadWeightsAndPoints.emplace_back(0.1294849661688697, -0.9491079123427585);
            quadWeightsAndPoints.emplace_back(0.2797053914892766, -0.7415311855993945);
            quadWeightsAndPoints.emplace_back(0.3818300505051189, -0.4058451513773972);
            quadWeightsAndPoints.emplace_back(0.4179591836734694, 0.0);
            quadWeightsAndPoints.emplace_back(0.3818300505051189, 0.4058451513773972);
            quadWeightsAndPoints.emplace_back(0.2797053914892766, 0.7415311855993945);
            quadWeightsAndPoints.emplace_back(0.1294849661688697, 0.9491079123427585);
        }

        default:{
            throw std::logic_error("Quadrature formula not available");
        }
    }
}

QuadPointData1D::QuadPointData1D(int n_quad_pts, int fe_order):n_quad_pts(n_quad_pts), fe_order(fe_order){
    GetWeightsAndPoints();
    basis_function reference(fe_order, -1.0, 1.0);
    quadratureWeightsAndGradients.reserve(n_quad_pts);
    for (auto it = quadWeightsAndPoints.begin(); it != quadWeightsAndPoints.end(); ++it)
    {
        std::vector<double> gradients_at_this_pt = reference.shape_function_gradients((*it).second);
        std::vector<double> functions_at_this_pt = reference.shape_function_values((*it).second);
        quadratureWeightsAndGradients.emplace_back((*it).first, gradients_at_this_pt);
        quadratureWeightsFunctionsAndGradients.emplace_back((*it).first, make_pair(functions_at_this_pt, gradients_at_this_pt));
    }
}

void Assemble::AssembleGlobalSystem(double kappa){
    unsigned int N_Cells = mesh.get_N_Cells();
    QuadPointData1D quadraturePointData(fe_order + 1, fe_order);
    vector<pair<double, vector<double>>>& quadratureWeightsAndGradients = quadraturePointData.quadratureWeightsAndGradients;
    
    for (unsigned int n_cell = 0; n_cell < N_Cells; n_cell++)
    {
        const cell_1d &cell = mesh.get_cell(n_cell);
        double x1 = cell.get_x1(), x2 = cell.get_x2();
        double mult = 2.0 / (x2 - x1);
        unsigned int N_Local_DOF = cell.get_N_DOF();
        double k_local[N_Local_DOF * N_Local_DOF];
        for (int i = 0; i < N_Local_DOF * N_Local_DOF; i++)
            k_local[i] = 0.0;

        for (unsigned int i_wq = 0; i_wq < quadraturePointData.n_quad_pts; i_wq++)
        {
            double quad_wt = quadratureWeightsAndGradients[i_wq].first;
            vector<double> &ref_grads = quadratureWeightsAndGradients[i_wq].second;
            double Nx[N_Local_DOF];

            assert((ref_grads.size() == N_Local_DOF));

            for (unsigned int i = 0; i < N_Local_DOF; i++)
                Nx[i] = ref_grads[i] * mult;

            for (unsigned int i = 0; i < N_Local_DOF; i++)
            {
                for (unsigned int j = 0; j < N_Local_DOF; j++)
                {
                    k_local[i * N_Local_DOF + j] += quad_wt * kappa * Nx[i] * Nx[j] / mult;
                }
            }
        }

        // for (int i = 0; i < N_Local_DOF; i++){
        //     for (int j = 0; j < N_Local_DOF; j++){
        //         cout << k_local[i * N_Local_DOF + j] << "\t";
        //     }
        //     cout << endl;
        // }
        // exit(0);
        
        unsigned int global_rows[N_Local_DOF];
        for (unsigned int i = 0; i < N_Local_DOF; i++){
            global_rows[i] = mesh.global_DOF(n_cell) + i;
        }

        matrix.add_local_matrix_to_global(k_local, global_rows);
    }
}

void Assemble::ApplyBoundaryCondition(int side, int cond){
    if(side == LEFT && cond == DIRICHLET){
        for (int i = matrix.RowPtr[0]; i < matrix.RowPtr[1]; i++){
            if(matrix.colIndex[i] == 0)
                matrix.entries[i] = 1;
            else
                matrix.entries[i] = 0;
        }
    }

    if(side == RIGHT && cond == DIRICHLET){
        for (int i = matrix.RowPtr[matrix.N_DOF - 1]; i < matrix.RowPtr[matrix.N_DOF]; i++){
            if(matrix.colIndex[i] == matrix.N_DOF - 1)
                matrix.entries[i] = 1;
            else
                matrix.entries[i] = 0;
        }
    }
}

void Assemble::ApplyBoundaryConditionSymmetric(int side, int cond, double val, vector<double>& rhs){
    int N_DOF = rhs.size();
    const int *RowPtr = matrix.RowPtr.data();
    const int *colIndex = matrix.colIndex.data();
    double *entries = matrix.entries.data();

    if (side == LEFT && cond == DIRICHLET)
    {
        for (int i = matrix.RowPtr[0]; i < matrix.RowPtr[1]; i++){
            if(matrix.colIndex[i] == 0)
                matrix.entries[i] = 1;
            else
                matrix.entries[i] = 0;
        }
        rhs[0] = val;
        for (int i = 0; i < N_DOF; i++){
            if(matrix.get(i, 0) == 0)
                break;
            for (int j = RowPtr[i]; j < RowPtr[i + 1]; j++)
            {
                if(colIndex[j] >= i)
                    break;
                if(colIndex[j] == 0){
                    rhs[i] -= entries[j] * val;
                    entries[j] = 0;
                }
            }
        }
    }

    if(side == RIGHT && cond == DIRICHLET){
        for (int i = matrix.RowPtr[matrix.N_DOF - 1]; i < matrix.RowPtr[matrix.N_DOF]; i++){
            if(matrix.colIndex[i] == matrix.N_DOF - 1)
                matrix.entries[i] = 1;
            else
                matrix.entries[i] = 0;
        }
        rhs[matrix.N_DOF - 1] = val;
        for (int i = 0; i < N_DOF; i++)
        {
            for (int j = RowPtr[i]; j < RowPtr[i + 1]; j++)
            {
                if (colIndex[j] == matrix.N_DOF - 1)
                {
                    if(i == colIndex[j])
                        break;
                    rhs[i] -= entries[j] * val;
                    entries[j] = 0;
                }
            }
        }
    }
}

void Assemble::AssembleGlobalSystemCD1D(double kappa, double convection_parameter){
    unsigned int N_Cells = mesh.get_N_Cells();
    QuadPointData1D quadraturePointData(fe_order + 1, fe_order);
    vector<pair<double, vector<double>>>& quadratureWeightsAndGradients = quadraturePointData.quadratureWeightsAndGradients;
    
    for (unsigned int n_cell = 0; n_cell < N_Cells; n_cell++)
    {
        const cell_1d &cell = mesh.get_cell(n_cell);
        double x1 = cell.get_x1(), x2 = cell.get_x2();
        double mult = 2.0 / (x2 - x1);
        unsigned int N_Local_DOF = cell.get_N_DOF();
        double k_local[N_Local_DOF * N_Local_DOF];
        for (int i = 0; i < N_Local_DOF * N_Local_DOF; i++)
            k_local[i] = 0.0;

        for (unsigned int i_wq = 0; i_wq < quadraturePointData.n_quad_pts; i_wq++)
        {
            double quad_wt = quadratureWeightsAndGradients[i_wq].first;
            vector<double> &ref_grads = quadratureWeightsAndGradients[i_wq].second;
            vector<double> &ref_funcs = quadraturePointData.quadratureWeightsFunctionsAndGradients[i_wq].second.first;
            double Nx[N_Local_DOF];

            assert((ref_grads.size() == N_Local_DOF));

            for (unsigned int i = 0; i < N_Local_DOF; i++)
                Nx[i] = ref_grads[i] * mult;

            for (unsigned int i = 0; i < N_Local_DOF; i++)
            {
                for (unsigned int j = 0; j < N_Local_DOF; j++)
                {
                    k_local[i * N_Local_DOF + j] += quad_wt * ((kappa * Nx[i] * Nx[j] / mult) - (convection_parameter * ref_funcs[i] * Nx[j]));
                }
            }
        }

        // for (int i = 0; i < N_Local_DOF; i++){
        //     for (int j = 0; j < N_Local_DOF; j++){
        //         cout << k_local[i * N_Local_DOF + j] << "\t";
        //     }
        //     cout << endl;
        // }
        // exit(0);
        
        unsigned int global_rows[N_Local_DOF];
        for (unsigned int i = 0; i < N_Local_DOF; i++){
            global_rows[i] = mesh.global_DOF(n_cell) + i;
        }

        matrix.add_local_matrix_to_global(k_local, global_rows);
    }
}

void Assemble::AssembleMassMatrix(){
    unsigned int N_Cells = mesh.get_N_Cells();
    QuadPointData1D quadraturePointData(fe_order + 1, fe_order);
    vector<pair<double, vector<double>>>& quadratureWeightsAndGradients = quadraturePointData.quadratureWeightsAndGradients;
    
    for (unsigned int n_cell = 0; n_cell < N_Cells; n_cell++)
    {
        const cell_1d &cell = mesh.get_cell(n_cell);
        double x1 = cell.get_x1(), x2 = cell.get_x2();
        double mult = 2.0 / (x2 - x1);
        unsigned int N_Local_DOF = cell.get_N_DOF();
        double k_local[N_Local_DOF * N_Local_DOF];
        for (int i = 0; i < N_Local_DOF * N_Local_DOF; i++)
            k_local[i] = 0.0;

        for (unsigned int i_wq = 0; i_wq < quadraturePointData.n_quad_pts; i_wq++)
        {
            double quad_wt = quadratureWeightsAndGradients[i_wq].first;
            vector<double> &ref_grads = quadratureWeightsAndGradients[i_wq].second;
            vector<double> &ref_funcs = quadraturePointData.quadratureWeightsFunctionsAndGradients[i_wq].second.first;
            double Nx[N_Local_DOF];

            assert((ref_grads.size() == N_Local_DOF));

            for (unsigned int i = 0; i < N_Local_DOF; i++)
                Nx[i] = ref_grads[i] * mult;

            for (unsigned int i = 0; i < N_Local_DOF; i++)
            {
                for (unsigned int j = 0; j < N_Local_DOF; j++)
                {
                    k_local[i * N_Local_DOF + j] += quad_wt * (ref_funcs[i] * ref_funcs[j]);
                }
            }
        }

        // for (int i = 0; i < N_Local_DOF; i++){
        //     for (int j = 0; j < N_Local_DOF; j++){
        //         cout << k_local[i * N_Local_DOF + j] << "\t";
        //     }
        //     cout << endl;
        // }
        // exit(0);
        
        unsigned int global_rows[N_Local_DOF];
        for (unsigned int i = 0; i < N_Local_DOF; i++){
            global_rows[i] = mesh.global_DOF(n_cell) + i;
        }

        mass_matrix.add_local_matrix_to_global(k_local, global_rows);
    }
}

void Assemble::GradFunctionForRK(const Eigen::VectorXd &u, double time, Eigen::VectorXd &du_dt){
    Eigen::VectorXd rhs(mesh.get_N_DOF());
    for (int i = 0; i < mesh.get_N_DOF(); i++)
    {
        rhs[i] = 0.0;
        for (int j = matrix.RowPtr[i]; j < matrix.RowPtr[i + 1]; j++)
            rhs[i] += matrix.entries[j] * u[matrix.colIndex[j]];

        rhs[i] *= -1.0;
    }

    Solver solver(100, 1e-6);
    solver.direct_solver(mass_matrix, rhs, du_dt, u.size());
}

// AssembleDG::AssembleDG(LegendreBasis &Lbasis, int order, const Mesh& mesh, double alpha) : 
// Lbasis(Lbasis), order(order), N_Local_DOF(order + 1), mesh(mesh), alpha(alpha), N_DOF((order+1)*mesh.get_N_Cells()) {

// }

// void AssembleDG::AssembleMkSLagrange(basis_function &lagrange){
//     const vector<double> &xi = lagrange.get_nodal_X();

//     M = Eigen::MatrixXd::Constant(N_Local_DOF, N_Local_DOF, 0.0);
//     S = Eigen::MatrixXd::Constant(N_Local_DOF, N_Local_DOF, 0.0);

//     QuadPointData1D qpdata(order + 1, order);

//     vector<pair<double, pair<vector<double>, vector<double>>>>& quadWtsFnsGrads = qpdata.quadratureWeightsFunctionsAndGradients;

//     for (int w_i = 0; w_i < qpdata.n_quad_pts; w_i++){
//         double quad_wt = quadWtsFnsGrads[w_i].first;
//         vector<double> &shape_fns = quadWtsFnsGrads[w_i].second.first;
//         vector<double> &shape_grads = quadWtsFnsGrads[w_i].second.second;

//         for (int i = 0; i < N_Local_DOF; i++)
//         {
//             for (int j = 0; j < N_Local_DOF; j++){
//                 M(i, j) += quad_wt * shape_fns[i] * shape_fns[j];
//                 S(i, j) += quad_wt * shape_fns[i] * shape_grads[j]; //j == 0 or j == N_Local_DOF - 1 ? fabs(j - (N_Local_DOF - 1)) : j
//             }
//         }
//     }
//     // for (int i = 0; i < N_Local_DOF; i++){
//     //     for (int j = 0; j < N_Local_DOF; j++){
//     //         cout << S(i, j) << "\t";
//     //     }
//     //     cout << endl;
//     // }
//     // exit(0);
//     // S = -1.0 * S;
// }

// void AssembleDG::AssembleMkS(){
//     vector<double> &r = Lbasis.LGLPoints;

//     Eigen::MatrixXd V(N_Local_DOF, N_Local_DOF), Vr(N_Local_DOF, N_Local_DOF);

//     for (int i = 0; i < N_Local_DOF; i++)
//     {
//         for (int j = 0; j < N_Local_DOF; j++){
//             V(i, j) = Lbasis.Pn(r[i], j, 0, 0);
//             Vr(i, j) = Lbasis.GradPn(r[i], j, 0, 0);
//         }
//     }

//     M = (V * V.transpose()).inverse();
//     Dr = Vr * V.inverse();
//     S = M * Dr;

//     hk.resize(mesh.get_N_Cells(), 0.0);

//     for (int i = 0; i < mesh.get_N_Cells(); i++){
//         const cell_1d &cell = mesh.get_cell(i);
//         hk[i] = (cell.get_x2() - cell.get_x1()) / 2.0;
//     }
// }

// void AssembleDG::AssembleRHS(Eigen::VectorXd& rhs, const Eigen::VectorXd& u){
//     int N_Cells = mesh.get_N_Cells();
//     double u_star;
//     for (int i = 0; i < N_Cells; i++)
//     {
//         if(!i){
//             u_star = ((u[N_Local_DOF - 1] + u[N_Local_DOF]) / 2) + ((1 - alpha) * (u[N_Local_DOF - 1] - u[N_Local_DOF]) / 2);
//             rhs[N_Local_DOF - 1] = 2 * M_PI * (u[N_Local_DOF - 1] - u_star);
//         }
//         else if(i == N_Cells - 1){
//             u_star = ((u[(N_Local_DOF * i) - 1] + u[N_Local_DOF * i]) / 2) +
//                      ((1 - alpha) * (u[(N_Local_DOF * i) - 1] - u[N_Local_DOF * i]) / 2);

//             rhs[N_Local_DOF * i] = -2 * M_PI * (u[N_Local_DOF * i] - u_star);
//         }
//         else{
//             u_star = ((u[(N_Local_DOF * i) - 1] + u[N_Local_DOF * i]) / 2) +
//                      ((1 - alpha) * (u[(N_Local_DOF * i) - 1] - u[N_Local_DOF * i]) / 2);

//             rhs[N_Local_DOF * i] = -2 * M_PI * (u[N_Local_DOF * i] - u_star);

//             u_star = ((u[N_Local_DOF * (i + 1) - 1] + u[N_Local_DOF * (i + 1)]) / 2) +
//                      ((1 - alpha) * (u[N_Local_DOF * (i + 1) - 1] - u[N_Local_DOF * (i + 1)]) / 2);

//             rhs[N_Local_DOF * (i + 1) - 1] = 2 * M_PI * (u[N_Local_DOF * (i + 1) - 1] - u_star);
//         }
//     }
//     // rhs[0] = rhs[N_Local_DOF * N_Cells - 1] = 0.0;
// }

// void AssembleDG::ApplyInlet(Eigen::VectorXd& rhs, const Eigen::VectorXd& u, double val){
//     double u_star = (val + u[0]) / 2 + ((1 - alpha) * (val - u[0]) / 2);

//     rhs[0] = -2 * M_PI * (u[0] - val); // val instead of u_star stopped divergence, dont know why!!
//     rhs[N_Local_DOF * mesh.get_N_Cells() - 1] = 0.0;
// }

// void AssembleDG::Getdu_dtForRK(const Eigen::VectorXd& u, const Eigen::VectorXd& rhs, Eigen::VectorXd &du_dt){
//     Eigen::MatrixXd M_inv = M.inverse();
//     Eigen::VectorXd uk(order + 1), rk(order + 1);

//     int N_Cells = mesh.get_N_Cells();
//     for (int i = 0; i < N_Cells; i++){
//         const cell_1d &cell = mesh.get_cell(i);
//         double hk = cell.get_x2() - cell.get_x1();

//         // memcpy(uk.data(), u.data() + (i * N_Local_DOF), N_Local_DOF * sizeof(double));
//         for (int j = 0; j < N_Local_DOF; j++)
//         {
//             uk[j] = u[i * N_Local_DOF + j];
//             rk[j] = rhs[i * N_Local_DOF + j];
//         }

//         Eigen::VectorXd dudt = (2 / hk) * ((M_inv * rk) - (2 * M_PI * (Dr * uk)));

//         for (int j = 0; j < N_Local_DOF; j++)
//             du_dt[i * N_Local_DOF + j] = dudt[j];
//     }
// }

// void AssembleDG::StepRK(Eigen::VectorXd &u, Eigen::VectorXd &rhs, double t, double dt, double val[]){
//     Eigen::VectorXd du_dt(N_DOF), k1(N_DOF), k2(N_DOF), k3(N_DOF), k4(N_DOF);
//     // rhs.resize(N_DOF);
//     // u.resize(N_DOF);
//     // k1.resize(N_DOF);
//     // k2.resize(N_DOF);
//     // k3.resize(N_DOF);
//     // k4.resize(N_DOF);
//     // du_dt.resize(N_DOF);

//     Eigen::VectorXd u_k1(N_DOF), u_k2(N_DOF), u_k3(N_DOF);

//     AssembleRHS(rhs, u);
//     ApplyInlet(rhs, u, val[0]);
//     Getdu_dtForRK(u, rhs, k1);
//     u_k1 = u + (0.5 * dt * k1);

//     AssembleRHS(rhs, u_k1);
//     ApplyInlet(rhs, u_k1, val[1]);
//     Getdu_dtForRK(u_k1, rhs, k2);
//     u_k2 = u + (0.5 * dt * k2);

//     AssembleRHS(rhs, u_k2);
//     ApplyInlet(rhs, u_k2, val[1]);
//     Getdu_dtForRK(u_k2, rhs, k3);
//     u_k3 = u + (dt * k3);

//     AssembleRHS(rhs, u_k3);
//     ApplyInlet(rhs, u_k3, val[2]);
//     Getdu_dtForRK(u_k2, rhs, k4);
//     u = (u + (1.0 / 6.0) * dt * (k1 + 2 * k2 + 2 * k3 + k4));
    
// }


// void AssembleDG::AssembleRHSWeak(Eigen::VectorXd& rhs, const Eigen::VectorXd& u){
//     int N_Cells = mesh.get_N_Cells();
//     double u_star;
//     for (int i = 0; i < N_Cells; i++)
//     {
//         if(!i){
//             u_star = ((u[N_Local_DOF - 1] + u[N_Local_DOF]) / 2) + ((1 - alpha) * (u[N_Local_DOF - 1] - u[N_Local_DOF]) / 2);
//             rhs[N_Local_DOF - 1] = -2.0 * M_PI * u_star;
//         }
//         else if(i == N_Cells - 1){
//             u_star = ((u[(N_Local_DOF * i) - 1] + u[N_Local_DOF * i]) / 2) +
//                      ((1 - alpha) * (u[(N_Local_DOF * i) - 1] - u[N_Local_DOF * i]) / 2);

//             rhs[N_Local_DOF * i] = 2.0 * M_PI * u_star;
//         }
//         else{
//             u_star = ((u[(N_Local_DOF * i) - 1] + u[N_Local_DOF * i]) / 2) +
//                      ((1 - alpha) * (u[(N_Local_DOF * i) - 1] - u[N_Local_DOF * i]) / 2);

//             rhs[N_Local_DOF * i] = 2.0 * M_PI * u_star;

//             u_star = ((u[N_Local_DOF * (i + 1) - 1] + u[N_Local_DOF * (i + 1)]) / 2) +
//                      ((1 - alpha) * (u[N_Local_DOF * (i + 1) - 1] - u[N_Local_DOF * (i + 1)]) / 2);

//             rhs[N_Local_DOF * (i + 1) - 1] = -2.0 * M_PI * u_star;
//         }
//     }
// }

// void AssembleDG::ApplyInletWeak(Eigen::VectorXd& rhs, const Eigen::VectorXd& u, double val){
//     double u_star = ((val + u[0]) / 2) + ((1 - alpha) * (val - u[0]) / 2);

//     rhs[0] = 2.0 * M_PI * u_star; // val instead of u_star stopped divergence, dont know why!!
//     rhs[N_Local_DOF * mesh.get_N_Cells() - 1] = -2.0 * M_PI * u[u.size() - 1];
// }

// void AssembleDG::Getdu_dtForRKWeak(const Eigen::VectorXd &u, const Eigen::VectorXd &rhs, Eigen::VectorXd &du_dt){
//     // cout << M << endl
//     //      << endl
//     //      << S << endl;
//     // exit(0);
//     Eigen::MatrixXd M_inv = M.inverse();
//     Eigen::MatrixXd M_invS = (M_inv * S.transpose());
//     Eigen::VectorXd uk(order + 1), rk(order + 1);

//     int N_Cells = mesh.get_N_Cells();
//     for (int i = 0; i < N_Cells; i++){
//         const cell_1d &cell = mesh.get_cell(i);
//         double hk = cell.get_x2() - cell.get_x1();

//         // memcpy(uk.data(), u.data() + (i * N_Local_DOF), N_Local_DOF * sizeof(double));
//         for (int j = 0; j < N_Local_DOF; j++)
//         {
//             uk[j] = u[i * N_Local_DOF + j];
//             rk[j] = rhs[i * N_Local_DOF + j];
//         }

//         Eigen::VectorXd dudt = (2 / hk) * ((M_inv * rk) + (2 * M_PI * (M_invS * uk)));

//         for (int j = 0; j < N_Local_DOF; j++)
//             du_dt[i * N_Local_DOF + j] = dudt[j];
//     }
// }

// void AssembleDG::StepRKWeak(Eigen::VectorXd &u, Eigen::VectorXd &rhs, double t, double dt, double val[]){
//     Eigen::VectorXd du_dt(N_DOF), k1(N_DOF), k2(N_DOF), k3(N_DOF), k4(N_DOF);
//     // rhs.resize(N_DOF);
//     // u.resize(N_DOF);
//     // k1.resize(N_DOF);
//     // k2.resize(N_DOF);
//     // k3.resize(N_DOF);
//     // k4.resize(N_DOF);
//     // du_dt.resize(N_DOF);

//     Eigen::VectorXd u_k1(N_DOF), u_k2(N_DOF), u_k3(N_DOF);

//     AssembleRHSWeak(rhs, u);
//     ApplyInletWeak(rhs, u, val[0]);
//     Getdu_dtForRKWeak(u, rhs, k1);
//     u_k1 = u + (0.5 * dt * k1);

//     AssembleRHSWeak(rhs, u_k1);
//     ApplyInletWeak(rhs, u_k1, val[1]);
//     Getdu_dtForRKWeak(u_k1, rhs, k2);
//     u_k2 = u + (0.5 * dt * k2);

//     AssembleRHSWeak(rhs, u_k2);
//     ApplyInletWeak(rhs, u_k2, val[1]);
//     Getdu_dtForRKWeak(u_k2, rhs, k3);
//     u_k3 = u + (dt * k3);

//     AssembleRHSWeak(rhs, u_k3);
//     ApplyInletWeak(rhs, u_k3, val[2]);
//     Getdu_dtForRKWeak(u_k2, rhs, k4);
//     u = (u + (1.0 / 6.0) * dt * (k1 + 2 * k2 + 2 * k3 + k4));
// }