#ifndef __POST_PROCESS_FOR_ASSIGNMENT__
#define __POST_PROCESS_FOR_ASSIGNMENT__
#include <vector>
#include <eigen3/Eigen/Dense>

void gnuplot(const std::vector<double> &x, const std::vector<double> &y);
void gnuplot(const std::vector<double> &x, const Eigen::VectorXd &y);
void write_vtk(const std::vector<double> &x, const std::vector<double> &y, const char *vtk_filename);
void write_vtk(const std::vector<double> &x, const Eigen::VectorXd &y, const char *filename);
void write_vtk_DG(const Mesh &mesh, const std::vector<double> &X, const Eigen::VectorXd &u, const char *filename);
void write_vtk_line(const std::vector<double> &x, const std::vector<double> &y, const char *filename);
void write_vtk_line(const std::vector<double> &x, const Eigen::VectorXd &y, const char *filename);
#endif