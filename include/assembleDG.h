#ifndef ASSEMBLE_FOR_DG_ASSIGNMENT
#define ASSEMBLE_FOR_DG_ASSIGNMENT
#include <mesh.h>
#include <basis_function.h>
#include <eigen3/Eigen/Dense>

class AssembleDG{
    const Mesh& mesh;
    int order, N_Local_DOF, N_DOF;
    double alpha;
    LegendreBasis& Lbasis;
    Eigen::MatrixXd M, S, Dr, S_tilde;
    std::vector<double> hk;

public:
    AssembleDG(int order, const Mesh& mesh, double alpha, LegendreBasis& Lbasis):order(order), mesh(mesh), alpha(alpha), Lbasis(Lbasis){
        
    }
    
    AssembleDG(LegendreBasis& Lbasis, int order, const Mesh& mesh, double alpha);

    void AssembleMkS();

    // void AssembleMkSLagrangeHeat(basis_function &, double);

    void AssembleMkSLagrange(basis_function &lagrange);

    void AssembleRHS(Eigen::VectorXd &rhs, const Eigen::VectorXd &u);

    void AssembleRHSWeak(Eigen::VectorXd &rhs, const Eigen::VectorXd &u, double, double rhs_ = 0.0);

    void ApplyInlet(Eigen::VectorXd &rhs, const Eigen::VectorXd &u, double val);

    void ApplyInletWeak(Eigen::VectorXd &rhs, const Eigen::VectorXd &u, double val, double coeff);

    void ApplyOutletWeak(Eigen::VectorXd &rhs, const Eigen::VectorXd &u, double val, double coeff);

    const Eigen::MatrixXd& getM() const{
        return M;
    }

    const Eigen::MatrixXd& getS() const{
        return S;
    }

    const Eigen::MatrixXd& getDr() const{
        return Dr;
    }

    const std::vector<double>& gethk() const{
        return hk;
    }

    void Getdu_dtForRK(const Eigen::VectorXd &u, const Eigen::VectorXd &rhs, Eigen::VectorXd &du_dt);

    void Getdu_dtForRKWeak(const Eigen::VectorXd &u, const Eigen::VectorXd &rhs, Eigen::VectorXd &du_dt);

    void StepRK(Eigen::VectorXd &u, Eigen::VectorXd &rhs, double t, double dt, double val[]);

    void StepRKWeak(Eigen::VectorXd &u, Eigen::VectorXd &rhs, double t, double dt, double val[]);
};
#endif