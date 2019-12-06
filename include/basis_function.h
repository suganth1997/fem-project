#ifndef __BASIS_FUNCTION_FOR_ASSIGNMENT__
#define __BASIS_FUNCTION_FOR_ASSIGNMENT__
#include <vector>

class basis_function
{
private:
    int ORDER;

    std::vector<double> X;

public:
    basis_function(int order, std::vector<double>& X):ORDER(order), X(X){

    }

    basis_function(int order, double x1, double x2);

    double evaluate(double x, const std::vector<double> &nodal_values);

    double gradient(double x);

    std::vector<double> shape_function_values(double x);

    std::vector<double> shape_function_gradients(double x);

    const std::vector<double> &get_nodal_X() const
    {
        return X;
    }

    int get_order() const{
        return ORDER;
    }
};

class LegendreBasis{
private:
    int ORDER;
public:
    std::vector<double> LGLPoints;

    LegendreBasis(int order);

    double a_val(int n, int alpha, int beta);

    double Pn(double r, int n, int alpha, int beta);

    double GradPn(double r, int n, int alpha, int beta);

};
#endif