#include <iostream>
#include <vector>
#include <map>
#include <mesh.h>
using namespace std;

cell_1d::cell_1d(double x_1, double x_2, unsigned int order, unsigned int id) : order(order){
    x1 = x_1;
    x2 = x_2;
    ID = id;
    N_DOF = order + 1;
    vertices.reserve(N_DOF);
    for (int i = 0; i < order + 1; i++){
        vertices.emplace_back(x1 + (x2 - x1) * i / order);
    }
    // for (double x = x_1; x < x_2 + ((x_2 - x_1) / (2 * order)); x += (x_2 - x_1) / order)
    //     vertices.emplace_back(x);
}

Mesh::Mesh(double l, double r, unsigned int n_cells, unsigned int order){
    left = l;
    right = r;
    N_Cells = n_cells;
    ORDER = order;
    N_DOF = N_Cells + 1 + (order - 1) * N_Cells;

    cells.reserve(N_Cells);
    local_to_global.reserve(N_Cells);
    double rln = (r - l) / N_Cells;
    for (int i = 0; i < N_Cells; i++){
        cells.emplace_back(l + i * rln, l + (i + 1) * rln, order, i);
        local_to_global.emplace_back(i == 0 ? i : local_to_global[i - 1] + order);
    }

    cells_shared_by_dof.reserve(N_DOF);
    for (int i = 0; i < N_DOF; i++)
        cells_shared_by_dof.emplace_back();

    for (int i = 0; i < N_Cells; i++){
        for (int j = 0; j < cells[i].get_N_DOF(); j++){
            unsigned int& N = cells_shared_by_dof[local_to_global[i] + j].N_Cells;
            cells_shared_by_dof[local_to_global[i] + j].cells[N] = i;
            N++;
        }
    }
}