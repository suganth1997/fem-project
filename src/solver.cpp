#include <iostream>
#include <vector>
#include <solver.h>
#include <cmath>
#include <ctime>
#include <cstring>
#include <sstream>
#include <umfpack_all.h>
#include <mkl.h>
#ifdef DEBUG
#define DEBUG_PRINT_VECTOR(x, n)             \
    for (int i = 0; i < n; i++) \
        cout << x[i] << "\t";   \
    cout << endl;
#define DEBUG_PRINT std::cout
#else
std::stringstream ss;
#define DEBUG_PRINT_VECTOR(x, n)
#define DEBUG_PRINT ss
#endif

using namespace std;

void Solver::solve_gauss_jacobi(Matrix &matrix, vector<double>& _rhs, vector<double>& _sol, int N_DOF, double omega = 1.0){
    if(omega < 1.0/10.0)
        throw logic_error("Omega value for jacobi is very bad");
    if (_sol.size() == 0)
        _sol.resize(N_DOF, 0.0);
    int *RowPtr = matrix.pattern.RowPtr.data();
    int *colIndex = matrix.pattern.colIndex.data();
    double *entries = matrix.entries.data();

    struct matrix_descr descrA;

    sparse_matrix_t csrA;    

    mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ZERO, N_DOF, N_DOF, RowPtr, RowPtr + 1, colIndex, entries);

    mkl_sparse_optimize(csrA);

    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

    int _max_iter = max_iter;
    double *rhs = _rhs.data();
    double *sol = _sol.data();
    vector<double> sol_prev, A_ii;
    sol_prev.resize(N_DOF, 0.0),
    A_ii.resize(N_DOF, 0.0);
    double buffer, residue = 1;
    for (int i = 0; i < N_DOF; i++)
    {
        A_ii[i] = 0.0;
        for (int j = RowPtr[i]; j < RowPtr[i + 1]; j++){
            if(i == colIndex[j]){
                A_ii[i] = entries[j];
                break;
            }
        }
    }

    int max_iter_ = max_iter;
    do
    {
        cblas_daxpby(N_DOF, 1.0, sol, 1, 0.0, sol_prev.data(), 1);
        mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrA, descrA, sol_prev.data(), 0.0, sol);

        for (int i = 0; i < N_DOF; i++)
        {
            buffer = sol[i];
            sol[i] = (rhs[i] - (sol[i] - A_ii[i] * sol_prev[i])) / A_ii[i];
            sol[i] = (1 - omega) * sol_prev[i] + omega * sol[i];
            sol_prev[i] = rhs[i] - buffer;

        }

        residue = cblas_dnrm2(N_DOF, sol_prev.data(), 1);

        DEBUG_PRINT << "Converged in " << max_iter_ - max_iter << " iterations" << endl
             << "Residual = " << residue << endl
             << endl;

    } while (residue > tol && max_iter--);

    cout << endl << "*******************************************************************" << endl;
    cout << "USING JACOBI SOLVER WITH OMEGA = " << omega << endl << endl;

    if (max_iter > 0)
    {
        cout << "Converged in " << max_iter_ - max_iter << " iterations" << endl
             << "Residual = " << residue << endl
             << endl;
    }
    else{
        cout << "Not converged" << endl;
    }

    cout << endl << "*******************************************************************" << endl;

}

void Solver::solve_gauss_siedel_mkl_naive(Matrix &matrix, vector<double>& _rhs, vector<double>& _sol, int N_DOF, double omega = 1.0){
    if(_sol.size() == 0)
        _sol.resize(N_DOF, 0.0);
    vector<int> &RowPtr = matrix.pattern.RowPtr;
    vector<int> &colIndex = matrix.pattern.colIndex;
    Matrix U(N_DOF, 0, matrix.pattern);
    for (int i = 0; i < N_DOF; i++){
        for (int j = RowPtr[i]; j < RowPtr[i + 1]; j++){
            if(colIndex[j] > i)
                U._get(i, colIndex[j]) = matrix.get(i, colIndex[j]);
        }
    }

    int *rowptr = RowPtr.data(), *colindex = colIndex.data();
    double *entries_U = U.entries.data();

    matrix_descr descrLD, descrU, descrA;
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
    descrU.type = SPARSE_MATRIX_TYPE_GENERAL;

    sparse_matrix_t csrLD, csrU, csrA;

    int *rowptrld, *colindld, *rowptrld_;

    sparse_index_base_t enum_ = SPARSE_INDEX_BASE_ZERO;

    double *entriesld;

    vector<double> Ux;
    Ux.resize(N_DOF, 0.0);
    // mkl_sparse_d_create_csr(&csrLD, SPARSE_INDEX_BASE_ZERO, N_DOF, N_DOF, rowptr, rowptr + 1, colindex, entries);
    mkl_sparse_d_create_csr(&csrU, SPARSE_INDEX_BASE_ZERO, N_DOF, N_DOF, rowptr, rowptr + 1, colindex, entries_U);
    mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ZERO, N_DOF, N_DOF, rowptr, rowptr + 1, colindex, matrix.entries.data());

    mkl_sparse_d_add(SPARSE_OPERATION_NON_TRANSPOSE, csrU, -1.0, csrA, &csrLD);

    descrLD.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
    descrLD.mode = SPARSE_FILL_MODE_LOWER;
    descrLD.diag = SPARSE_DIAG_NON_UNIT;
    int max_iter_ = max_iter;
    double residue = 1;

    do
    {
        vector<double> sol_old = _sol;
        mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, csrU, descrU, _sol.data(), 0.0, Ux.data());
        cblas_daxpy(N_DOF, 1.0, _rhs.data(), 1, Ux.data(), 1);
        
        mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrLD, descrLD, Ux.data(), _sol.data());
        DEBUG_PRINT_VECTOR(_sol, N_DOF);
        // exit(0);
        // mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrA, descrA, _sol.data(), 0.0, Ux.data());

        cblas_daxpy(N_DOF, -1.0, _sol.data(), 1, sol_old.data(), 1);
        // cblas_daxpy(N_DOF, -1.0, _rhs.data(), 1, Ux.data(), 1);

        residue = cblas_dnrm2(N_DOF, sol_old.data(), 1);
        DEBUG_PRINT << "Converged in " << max_iter_ - max_iter << " iterations" << endl
             << "Residual = " << residue << endl
             << endl;
    } while (residue > tol && max_iter--);

    

    cout << endl << "*******************************************************************" << endl;
    cout << "USING GAUSS SIEDEL SOLVER WITH OMEGA = " << omega << endl << endl;

    if (max_iter > 0)
    {
        cout << "Converged in " << max_iter_ - max_iter << " iterations" << endl
             << "Residual = " << residue << endl
             << endl;
    }
    else{
        cout << "Not converged" << endl;
    }
    cout << endl
         << "*******************************************************************" << endl;
    // mkl_sparse_d_export_csr(csrLD, &enum_, &N_DOF, &N_DOF, &rowptrld, &rowptrld_, &colindld, &entriesld);

    // double LD[N_DOF][N_DOF];
    // for(int i=0; i<N_DOF; i++){
    //     for(int j = 0; j<N_DOF; j++){
    //         LD[i][j] = 0.0;
    //     }

    //     for(int j=rowptrld[i]; j<rowptrld_[i]; j++){
    //         LD[i][colindld[j]] = entriesld[j];
    //     }

    //     for(int j=0; j<N_DOF; j++){
    //         cout << LD[i][j] << "\t";
    //     }
    //     cout << endl;
    // }

    // exit(0);
}

void Solver::solve_gauss_siedel(Matrix &matrix, std::vector<double> &rhs, std::vector<double> &sol, int N_DOF, double omega = 1.0){
    if(sol.size() == 0)
        sol.resize(N_DOF, 0.0);

    // vector<int> &RowPtr_ = matrix.pattern.RowPtr;
    // vector<int> &colIndex_ = matrix.pattern.colIndex;

    // int *rowptr = RowPtr_.data(), *colindex = colIndex_.data();
    // double *entries_ = matrix.entries.data();

    matrix_descr descrA;
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

    sparse_matrix_t csrA;

    vector<double> x0 = sol;
    int max_iter_ = max_iter;
    double norm = 0.0;
    const int *RowPtr = matrix.RowPtr.data();
    const int *colIndex = matrix.colIndex.data();
    const double *entries = matrix.entries.data();
    do
    {
        norm = 0;
        for( int i = 0 ; i < N_DOF; i++)
        {   
            x0[i] = sol[i];
            sol[i] = rhs[i];
            double A_ii = 0.0;
            for (int j = RowPtr[i]; j < RowPtr[i + 1]; j++)
                if(i != colIndex[j])
                    sol[i] -= entries[j]*sol[colIndex[j]];
                else
                {
                    A_ii = entries[j];
                }

            sol[i] = (1 - omega) * x0[i] + omega * sol[i] / A_ii;
            norm += (x0[i] - sol[i]) * (x0[i] - sol[i]);
        }

        DEBUG_PRINT << "Converged in " << max_iter_ - max_iter << " iterations" << endl
             << "Residual = " << norm << endl
             << endl;
        // PRINT(sol, N_DOF);
        // exit(0);
        norm = sqrt(norm);
        //residual = residual_norm(matrix,b,x);
    }

    while (norm > tol && max_iter--);

    cout << endl << "*******************************************************************" << endl;
    cout << "USING GAUSS SIEDEL SOLVER WITH OMEGA = " << omega << endl << endl;

    if (max_iter > 0)
    {
        cout << "Converged in " << max_iter_ - max_iter << " iterations" << endl
             << "Residual = " << norm << endl
             << endl;
    }
    else{
        cout << "Not converged" << endl;
    }

    cout << endl << "*******************************************************************" << endl;
}

void Solver::direct_solver(Matrix &matrix, std::vector<double> &rhs, std::vector<double> &sol, unsigned int N_DOF){
    if(sol.size() != N_DOF)
        sol.resize(N_DOF, 0.0);
    
    const int *RowPtr, *colIndex;
    double *entries;

    RowPtr = matrix.RowPtr.data();
    colIndex = matrix.colIndex.data();
    entries = matrix.entries.data();

    void *Symbolic, *Numeric;

    int ret_flag = umfpack_di_symbolic(N_DOF, N_DOF, RowPtr, colIndex, entries, &Symbolic, nullptr, nullptr);

    if(ret_flag != 0){
        cout << "Error in umfpack_di_symbolic" << endl;
        exit(0);
    }
    UMFPACK_return(ret_flag);
    ret_flag = umfpack_di_numeric(RowPtr, colIndex, entries, Symbolic, &Numeric, nullptr, nullptr);
    umfpack_di_free_symbolic(&Symbolic);

    if(ret_flag != 0){
        cout << "Error in umfpack_di_numeric" << endl;
        exit(0);
    }
    UMFPACK_return(ret_flag);
    umfpack_di_solve(UMFPACK_At, RowPtr, colIndex, entries, sol.data(), rhs.data(), Numeric, nullptr, nullptr);
    umfpack_di_free_numeric(&Numeric);

    if(ret_flag != 0){
        cout << "Error in umfpack_di_solve" << endl;
        exit(0);
    }

    UMFPACK_return(ret_flag);
    cout << endl
         << "DIRECT SOLVER SUCCESSFUL" << endl;
}

void Solver::direct_solver(Matrix &matrix, Eigen::VectorXd &rhs, Eigen::VectorXd &sol, unsigned int N_DOF){
    const int *RowPtr, *colIndex;
    double *entries;

    RowPtr = matrix.RowPtr.data();
    colIndex = matrix.colIndex.data();
    entries = matrix.entries.data();

    void *Symbolic, *Numeric;

    int ret_flag = umfpack_di_symbolic(N_DOF, N_DOF, RowPtr, colIndex, entries, &Symbolic, nullptr, nullptr);

    if(ret_flag != 0){
        cout << "Error in umfpack_di_symbolic" << endl;
        exit(0);
    }
    UMFPACK_return(ret_flag);
    ret_flag = umfpack_di_numeric(RowPtr, colIndex, entries, Symbolic, &Numeric, nullptr, nullptr);
    umfpack_di_free_symbolic(&Symbolic);

    if(ret_flag != 0){
        cout << "Error in umfpack_di_numeric" << endl;
        exit(0);
    }
    UMFPACK_return(ret_flag);
    umfpack_di_solve(UMFPACK_At, RowPtr, colIndex, entries, sol.data(), rhs.data(), Numeric, nullptr, nullptr);
    umfpack_di_free_numeric(&Numeric);

    if(ret_flag != 0){
        cout << "Error in umfpack_di_solve" << endl;
        exit(0);
    }

    UMFPACK_return(ret_flag);
    cout << endl
         << "DIRECT SOLVER SUCCESSFUL" << endl;
}

void Solver::solve_cg(Matrix &matrix, std::vector<double> &rhs, std::vector<double> &sol, int N_DOF){
    if(sol.size() != N_DOF)
        sol.resize(N_DOF, 0);

    vector<double> _b = rhs;

    int *RowPtr = matrix.pattern.RowPtr.data();
    int *colIndex = matrix.pattern.colIndex.data();
    double *entries = matrix.entries.data();

    struct matrix_descr descrA;

    sparse_matrix_t csrA;    

    mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ZERO, N_DOF, N_DOF, RowPtr, RowPtr + 1, colIndex, entries);

    mkl_sparse_optimize(csrA);

    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

    vector<double> _p, _Ap;
    _p.resize(N_DOF, 0);
    _Ap.resize(N_DOF, 0);
    double *r = _b.data(), *x = sol.data(), *p = _p.data(), *Ap = _Ap.data(); /////////////////////// Ap is residual

    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, csrA, descrA, x, 1.0, r);

    double res_norm = cblas_dnrm2(N_DOF, r, 1);
    int max_iter_ = max_iter = N_DOF;

    if(res_norm < tol)
        return;

    memcpy(p, r, N_DOF * sizeof(double));

    double alpha = 1.0;
    while (res_norm > tol && max_iter--)
    {
        double rTr = cblas_ddot(N_DOF, r, 1, r, 1);

        mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrA, descrA, p, 0.0, Ap);

        double pTAp = cblas_ddot(N_DOF, p, 1, Ap, 1);

        alpha = rTr / pTAp;

        cblas_daxpy(N_DOF, alpha, p, 1, x, 1);

        cblas_daxpy(N_DOF, -1.0 * alpha, Ap, 1, r, 1);

        res_norm = cblas_dnrm2(N_DOF, r, 1);

        if (res_norm < tol)
            break;

        double r_T_r = cblas_ddot(N_DOF, r, 1, r, 1);

        double beta = r_T_r / rTr;

        cblas_daxpby(N_DOF, 1.0, r, 1, beta, p, 1);

        DEBUG_PRINT << "Residual Norm = " << res_norm << " at iteration " << max_iter_ - max_iter << endl;
    }

    cout << endl << "*******************************************************************" << endl;
    cout << "USING CONJUGATE GRADIENT SOLVER" << endl << endl;

    if (max_iter > 0)
    {
        cout << "Converged in " << max_iter_ - max_iter << " iterations";
    }
    else{
        cout << "Not converged" << endl;
    }

    cout << endl << "*******************************************************************" << endl;
}