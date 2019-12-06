#include <iostream>
#include <string>
#include <algorithm>
#include <iomanip>
#include <matrix.h>
using namespace std;
const int maximum_dof_until_which_debug_print_function_should_work = 25;

void sparsity_pattern::compute_sparsity(){
    unsigned int N_DOF = mesh.get_N_DOF();
    RowPtr.reserve(mesh.get_N_DOF() + 1);
    RowPtr.emplace_back(0);

    nnz = 0;
    for (int i = 0; i < N_DOF; i++){
        const DOF_Info &dof_info = mesh.get_cell_sharing_info(i);
        unsigned int n_this_row = dof_info.N_Cells == 1 ? N_BaseFunctions : N_BaseFunctions * dof_info.N_Cells - 1;
        RowPtr.emplace_back(RowPtr[i] + n_this_row);
        nnz += n_this_row;
    }

    colIndex.reserve(nnz);

    for (int i = 0; i < N_DOF; i++){
        const DOF_Info &dof_info = mesh.get_cell_sharing_info(i);
        if(dof_info.N_Cells == 1){
            for (int j = 0; j < N_BaseFunctions; j++)
                colIndex.emplace_back(mesh.global_DOF(dof_info.cells[0]) + j);
        }
        else if(dof_info.N_Cells == 2){
            for (int j = 0; j < 2 * N_BaseFunctions - 1; j++)
                colIndex.emplace_back(mesh.global_DOF(dof_info.cells[0]) + j);
        }
    }
}


void sparsity_pattern::print_sparsity_pattern() const{
    unsigned int N_DOF = mesh.get_N_DOF();
    if(N_DOF > maximum_dof_until_which_debug_print_function_should_work)
        throw logic_error("Not feasible to debug print with huge DOFs");
    
    char A[N_DOF][N_DOF];
    for (int i = 0; i < N_DOF; i++)
        for (int j = 0; j < N_DOF; j++)
            A[i][j] = '-';

    for (int i = 0; i < N_DOF; i++)
    {
        for (int j = RowPtr[i]; j < RowPtr[i + 1]; j++)
        {
            A[i][colIndex[j]] = '+';
        }
    }

    for (int i = 0; i < N_DOF; i++){
        for (int j = 0; j < N_DOF; j++){
            cout << A[i][j] << "\t";
        }
        cout << endl;
    }
}

void Matrix::add_local_matrix_to_global(double k_local[], unsigned int global_rows[]){

    for (unsigned int i = 0; i < N_Local; i++)
    {
        unsigned int local_j = 0;
        for (unsigned int j = RowPtr[global_rows[i]]; j < RowPtr[global_rows[i] + 1]; j++)
        {
            if(colIndex[j] == global_rows[local_j])
            {
                entries[j] += k_local[i * N_Local + local_j];
                local_j++;
            }
        }
    }
}

double Matrix::get(unsigned int _i, unsigned int j){
    for (int i = RowPtr[_i]; i < RowPtr[_i + 1]; i++){
        if(colIndex[i] == j)
            return entries[i];
    }
    return 0.0;
}

double& Matrix::_get(unsigned int _i, unsigned int j){
    for (int i = RowPtr[_i]; i < RowPtr[_i + 1]; i++){
        if(colIndex[i] == j)
            return entries[i];
    }
    cout << "Error at" << endl
         << "i = " << _i << "j = " << j << endl;
    throw logic_error("You are accesing zero entries so I can't return a reference");
}

void Matrix::print(string key) const{
    if(N_DOF > maximum_dof_until_which_debug_print_function_should_work){
        cout << "Not feasible to print huge matrices" << endl;
        return;
    }

    transform(key.begin(), key.end(), key.begin(), ::tolower);
    double A[N_DOF][N_DOF];

    for (int i = 0; i < N_DOF; i++)
        for (int j = 0; j < N_DOF; j++)
            A[i][j] = 0;

    for (int i = 0; i < N_DOF; i++)
    {
        for (int j = RowPtr[i]; j < RowPtr[i + 1]; j++)
        {
            A[i][colIndex[j]] = entries[j];
        }
    }
    cout << setprecision(3);
    if(!key.compare("python")){
        cout << "np.array([";
        for (int i = 0; i < N_DOF; i++)
        {
            cout << "[";
            for (int j = 0; j < N_DOF; j++){
                cout << A[i][j] << ", ";
            }
            cout << "]," << endl;
        }
        cout << "])" << endl;
    }
    else{
        for (int i = 0; i < N_DOF; i++)
        {
            for (int j = 0; j < N_DOF; j++)
                cout << A[i][j] << "\t";
            if(!key.compare("matlab"))
                cout << ";";
            cout << endl;
        }
    }
}