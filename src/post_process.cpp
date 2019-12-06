#include <sstream>
#include <string>
#include <fstream>
#include <mesh.h>
#include <post_process.h>
#include <eigen3/Eigen/Dense>

#ifndef GET_VARIABLE_NAME
#define GET_VARIABLE_NAME(Variable) (#Variable)
#endif

using namespace std;
void gnuplot(const vector<double> &x, const vector<double> &y)
{
    const char p[] = "_GNUPLOT_TEMP_FILE";
    ofstream file(p);
    for (int i = 0; i < x.size(); i++)
        file << x[i] << "\t" << y[i] << endl;

    system("gnuplot -p -e 'plot \"_GNUPLOT_TEMP_FILE\" using 1:2 with lines'");
    system("rm _GNUPLOT_TEMP_FILE");
}

void gnuplot(const vector<double> &x, const Eigen::VectorXd &y)
{
    const char p[] = "_GNUPLOT_TEMP_FILE";
    ofstream file(p);
    for (int i = 0; i < x.size(); i++)
        file << x[i] << "\t" << y(i) << endl;

    system("gnuplot -p -e 'plot \"_GNUPLOT_TEMP_FILE\" using 1:2 with lines'");
    system("rm _GNUPLOT_TEMP_FILE");
}

void write_vtk(const std::vector<double> &x, const std::vector<double> &y, const char *filename){
    stringstream ss;
    ss << "VTK/" << filename << ".vtk";
    system("mkdir -p VTK");

    ofstream file(ss.str().c_str());

    file << "# vtk DataFile Version 3.0" << endl
         << "vtk output" << endl
         << "ASCII" << endl
         << "DATASET UNSTRUCTURED_GRID" << endl
         << endl;

    file << "POINTS " << x.size() << " float" << endl;

    for (int i = 0; i < x.size(); i++){
        file << x[i] << " 0 0" << endl;
    }

    file << endl;

    file << "CELLS " << x.size() / 2 << " " << (x.size() / 2) * 3 << endl;

    int cell_i = 0;

    for (int i = 0; i < x.size() / 2; i++){
        file << "2 " << cell_i << " " << cell_i + 1 << endl;
        cell_i += 2;
    }

    file << "CELL_TYPES " << x.size() / 2 << endl;

    for (int i = 0; i < x.size() / 2; i++)
        file << "3 ";

    file << endl;

    file << "POINT_DATA " << x.size() << endl;
    file << "SCALARS u float" << endl
         << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < y.size(); i++)
        file << y[i] << " ";

    file << endl;
}

void write_vtk(const std::vector<double> &x, const Eigen::VectorXd &y, const char *filename){
    stringstream ss;
    ss << "VTK/" << filename << ".vtk";
    system("mkdir -p VTK");

    ofstream file(ss.str().c_str());

    file << "# vtk DataFile Version 3.0" << endl
         << "vtk output" << endl
         << "ASCII" << endl
         << "DATASET UNSTRUCTURED_GRID" << endl
         << endl;

    file << "POINTS " << x.size() << " float" << endl;

    for (int i = 0; i < x.size(); i++){
        file << x[i] << " 0 0" << endl;
    }

    file << endl;

    file << "CELLS " << x.size() / 2 << " " << (x.size() / 2) * 3 << endl;

    int cell_i = 0;

    for (int i = 0; i < x.size() / 2; i++){
        file << "2 " << cell_i << " " << cell_i + 1 << endl;
        cell_i += 2;
    }

    file << "CELL_TYPES " << x.size() / 2 << endl;

    for (int i = 0; i < x.size() / 2; i++)
        file << "3 ";

    file << endl;

    file << "POINT_DATA " << x.size() << endl;
    file << "SCALARS u float" << endl
         << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < x.size(); i++)
        file << y[i] << " ";

    file << endl;
}

void write_vtk_DG(const Mesh& mesh, const vector<double>& X, const Eigen::VectorXd& u, const char *filename){
    vector<double> Y;
    int N_Local_DOF = mesh.get_order();
    Y.reserve(mesh.get_N_Cells() * N_Local_DOF);
    for (int i = 0; i < mesh.get_N_Cells(); i++)
    {
        Y.push_back(u(i * N_Local_DOF));
        Y.push_back(u((i + 1) * N_Local_DOF - 1));
    }
    write_vtk(X, Y, filename);
}

void write_vtk_line(const std::vector<double> &x, const std::vector<double> &y, const char *filename){
    stringstream ss;
    ss << "VTK/" << filename << ".vtk";
    system("mkdir -p VTK");

    ofstream file(ss.str().c_str());

    file << "# vtk DataFile Version 3.0" << endl
         << "vtk output" << endl
         << "ASCII" << endl
         << "DATASET UNSTRUCTURED_GRID" << endl
         << endl;

    file << "POINTS " << x.size() << " float" << endl;

    for (int i = 0; i < x.size(); i++){
        file << x[i] << " 0 0" << endl;
    }

    file << endl;

    file << "CELLS " << x.size() - 1 << " " << (x.size() - 1) * 3 << endl;


    for (int i = 0; i < x.size() - 1; i++){
        file << "2 " << i << " " << i + 1 << endl;
    }

    file << "CELL_TYPES " << x.size() / 2 << endl;

    for (int i = 0; i < x.size() - 1; i++)
        file << "3 ";

    file << endl;

    file << "POINT_DATA " << x.size() << endl;
    file << "SCALARS u float" << endl
         << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < y.size(); i++)
        file << y[i] << " ";

    file << endl;
}


void write_vtk_line(const std::vector<double> &x, const Eigen::VectorXd &y, const char *filename){
    stringstream ss;
    ss << "VTK/" << filename << ".vtk";
    system("mkdir -p VTK");

    ofstream file(ss.str().c_str());

    file << "# vtk DataFile Version 3.0" << endl
         << "vtk output" << endl
         << "ASCII" << endl
         << "DATASET UNSTRUCTURED_GRID" << endl
         << endl;

    file << "POINTS " << x.size() << " float" << endl;

    for (int i = 0; i < x.size(); i++){
        file << x[i] << " 0 0" << endl;
    }

    file << endl;

    file << "CELLS " << x.size() - 1 << " " << (x.size() - 1) * 3 << endl;


    for (int i = 0; i < x.size() - 1; i++){
        file << "2 " << i << " " << i + 1 << endl;
    }

    file << "CELL_TYPES " << x.size() - 1 << endl;

    for (int i = 0; i < x.size() - 1; i++)
        file << "3 ";

    file << endl;

    file << "POINT_DATA " << x.size() << endl;
    file << "SCALARS u float" << endl
         << "LOOKUP_TABLE default" << endl;

    for (int i = 0; i < y.size(); i++)
        file << y[i] << " ";

    file << endl;
}