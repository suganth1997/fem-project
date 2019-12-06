#include <vector>
#include <matrix.h>
#include <eigen3/Eigen/Dense>

class Solver{
    int max_iter;
    double tol;

public:
    Solver(int max_iter, double tol):max_iter(max_iter), tol(tol){

    }

    void solve_gauss_jacobi(Matrix &matrix, std::vector<double> &rhs, std::vector<double> &sol, int N_DOF, double omega);

    void solve_gauss_siedel_mkl_naive(Matrix &matrix, std::vector<double> &rhs, std::vector<double> &sol, int N_DOF, double omega);

    void solve_gauss_siedel(Matrix &matrix, std::vector<double> &rhs, std::vector<double> &sol, int N_DOF, double omega);

    void solve_cg(Matrix &matrix, std::vector<double> &rhs, std::vector<double> &sol, int N_DOF);

    void direct_solver(Matrix &Matrix, std::vector<double> &rhs, std::vector<double> &sol, unsigned int N_DOF);

    void direct_solver(Matrix &Matrix, Eigen::VectorXd &rhs, Eigen::VectorXd &sol, unsigned int N_DOF);

};

