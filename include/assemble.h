#ifndef __ASSEMBLE_FOR_ASSIGNMENT__
#define __ASSEMBLE_FOR_ASSIGNMENT__

#include <matrix.h>
#include <mesh.h>
#include <basis_function.h>
#include <vector>
#include <eigen3/Eigen/Dense>

enum BoundPoint
{
    LEFT,
    RIGHT
};

enum BoundaryCondition
{
    DIRICHLET,
    NEUMANN
};

class Assemble{
    int fe_order;
    const Mesh &mesh;
    Matrix &matrix, &mass_matrix;

public:
    Assemble(const Mesh &mesh, Matrix &matrix, Matrix &mass_matrix, int fe_order) : mesh(mesh), matrix(matrix), mass_matrix(mass_matrix), fe_order(fe_order){

    }

    void AssembleGlobalSystem(double);

    void AssembleGlobalSystemCD1D(double, double);

    void AssembleMassMatrix();

    void GradFunctionForRK(const Eigen::VectorXd &u, double time, Eigen::VectorXd &du_dt);

    void ApplyBoundaryCondition(int side, int cond);

    void ApplyBoundaryConditionSymmetric(int side, int cond, double val, std::vector<double> &rhs);
};

class QuadPointData1D{
    int n_quad_pts, fe_order;
    std::vector<std::pair<double, std::vector<double>>> quadratureWeightsAndGradients;
    std::vector<std::pair<double, std::pair<std::vector<double>, std::vector<double>>>> quadratureWeightsFunctionsAndGradients;
    std::vector<std::pair<double, double>> quadWeightsAndPoints;

    void GetWeightsAndPoints();

public:
    QuadPointData1D(int gaussian_order, int fe_order);

    friend class Assemble;
    friend class AssembleDG;
};

// class AssembleDG{
//     const Mesh& mesh;
//     int order, N_Local_DOF, N_DOF;
//     double alpha;
//     LegendreBasis &Lbasis;
//     Eigen::MatrixXd M, S, Dr;
//     std::vector<double> hk;

// public:
//     AssembleDG(LegendreBasis& Lbasis, int order, const Mesh& mesh, double alpha);

//     void AssembleMkS();

//     void AssembleMkSLagrange(basis_function &lagrange);

//     void AssembleRHS(Eigen::VectorXd &rhs, const Eigen::VectorXd &u);

//     void AssembleRHSWeak(Eigen::VectorXd &rhs, const Eigen::VectorXd &u);

//     void ApplyInlet(Eigen::VectorXd &rhs, const Eigen::VectorXd &u, double val);

//     void ApplyInletWeak(Eigen::VectorXd &rhs, const Eigen::VectorXd &u, double val);

//     const Eigen::MatrixXd& getM() const{
//         return M;
//     }

//     const Eigen::MatrixXd& getS() const{
//         return S;
//     }

//     const Eigen::MatrixXd& getDr() const{
//         return Dr;
//     }

//     const std::vector<double>& gethk() const{
//         return hk;
//     }

//     void Getdu_dtForRK(const Eigen::VectorXd &u, const Eigen::VectorXd &rhs, Eigen::VectorXd &du_dt);

//     void Getdu_dtForRKWeak(const Eigen::VectorXd &u, const Eigen::VectorXd &rhs, Eigen::VectorXd &du_dt);

//     void StepRK(Eigen::VectorXd &u, Eigen::VectorXd &rhs, double t, double dt, double val[]);

//     void StepRKWeak(Eigen::VectorXd &u, Eigen::VectorXd &rhs, double t, double dt, double val[]);
// };
#endif