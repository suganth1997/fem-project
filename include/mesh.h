#ifndef __MESHING_FOR_ASSIGNMENT__
#define __MESHING_FOR_ASSIGNMENT__
#include <vector>
#include <stdexcept>

struct DOF_Info{
    unsigned int N_Cells = 0, cells[2] = {0};
    // DOF_Info(unsigned int N_Cells, unsigned int cell_1, unsigned int cell_2 = 0) : N_Cells(N_Cells){
    //     cell_ids[0] = cell_1;
    //     cell_ids[1] = cell_2;
    // }
};

class cell_1d
{
private:
    unsigned int ID, order, N_DOF;
    double x1, x2;
    std::vector<double> vertices;

public:
    cell_1d(double x_1, double x_2, unsigned int order, unsigned int id);

    double get_x1() const{
        return x1;
    }

    double get_x2() const{
        return x2;
    }

    double get_N_DOF() const{
        return N_DOF;
    }

    const std::vector<double>& get_vertices() const{
        return vertices;
    }
};


class Mesh
{
private:
    unsigned int ORDER, N_Cells, N_DOF;
    double left, right;
    std::vector<cell_1d> cells;
    std::vector<unsigned int> local_to_global; // map[cell_id] = beginindex
    std::vector<DOF_Info> cells_shared_by_dof;

public:
    Mesh(double l, double r, unsigned int N_Cells, unsigned int order);

    unsigned int get_N_Cells() const{
        return N_Cells;
    }

    unsigned int get_N_DOF() const{
        return N_DOF;
    }

    unsigned int global_DOF(int i) const{
        return local_to_global[i];
    }

    unsigned int get_order() const{
        return ORDER;
    }

    const cell_1d &get_cell(int i) const{
        if(i < 0 || i >= N_Cells)
            throw std::logic_error("Invalid Index");
        
        return cells[i];
    }

    const DOF_Info& get_cell_sharing_info(unsigned int dof) const{
        return cells_shared_by_dof[dof];
    }
};
#endif