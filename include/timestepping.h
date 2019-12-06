#include <functional>
#include <vector>
#include <assemble.h>
#include <eigen3/Eigen/Dense>

using GradientForRK = std::function<void (const Eigen::VectorXd &u, double time, Eigen::VectorXd &du_dt)>;

// typedef void(Assemble::*)(const Eigen::VectorXd, double, Eigen::VectorXd) GradientForRK;

void StepRK4(Eigen::VectorXd &u_init, double curr_time, double dt, int N_DOF, Assemble& assemble);