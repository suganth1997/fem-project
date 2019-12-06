#ifndef __ASSIGNMENT_MATRIX___
#define __ASSIGNMENT_MATRIX___
#include <vector>
#include <mesh.h>

class sparsity_pattern{
    unsigned int nnz, N_BaseFunctions;
    std::vector<int> RowPtr;
    std::vector<int> colIndex;
    const Mesh &mesh;

    std::vector<int>& get_RowPtr(){
        return RowPtr;
    }

    std::vector<int>& get_colIndex(){
        return colIndex;
    }

public:
    sparsity_pattern(const Mesh &mesh, unsigned int N) : mesh(mesh), N_BaseFunctions(N){

    }

    unsigned int get_nnz() const{
        return nnz;
    }

    void compute_sparsity();

    void compute_sparsityDifferent();

    void print_sparsity_pattern() const;

    const std::vector<int>& get_RowPtr() const{
        return RowPtr;
    }

    const std::vector<int>& get_colIndex() const{
        return colIndex;
    }

    friend class Matrix;
    friend class Solver;
};

class Matrix
{
public:
    unsigned int N_DOF, N_Local;
    std::vector<double> entries;
    const std::vector<int> RowPtr, colIndex;
    sparsity_pattern &pattern;

    void add_local_matrix_to_global(double k_local[], unsigned int global_rows[]);

    void clear(){
        for (int i = 0; i < entries.size(); i++)
            entries[i] = 0.0;
    }


    Matrix(unsigned int N_DOF, unsigned int N_Local, sparsity_pattern& pattern)
    :N_DOF(N_DOF), N_Local(N_Local), pattern(pattern), RowPtr(pattern.get_RowPtr()), colIndex(pattern.get_colIndex()){
        entries.resize(pattern.get_nnz(), 0);
    }

    Matrix(unsigned int N_DOF, const std::vector<int> &RowPtr, const std::vector<int> &colIndex, const std::vector<double> &entries, sparsity_pattern &pattern)
    :N_DOF(N_DOF), RowPtr(RowPtr), colIndex(colIndex), entries(entries), pattern(pattern){

    }

    void print(std::string key = "") const;

    double get(unsigned int i, unsigned int j);

    double& _get(unsigned int i, unsigned int j);

    const std::vector<int>& get_RowPtr() const{
        return RowPtr;
    }

    const std::vector<int>& get_colIndex() const{
        return colIndex;
    }

    const std::vector<double>& get_entries(){
        return entries;
    }

    friend class Assemble;
    friend class Solver;
};
#endif