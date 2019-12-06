#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <sstream>
#include <exception>
#include <fstream>
#include <basis_function.h>
using namespace std;

basis_function::basis_function(int order, double x1, double x2){
    if(x1 == x2)
        throw logic_error("Cannot create basis function");

    if(x2 < x1)
        swap(x1, x2);

    ORDER = order;
    X.push_back(x1);
    for (int i = 0; i < order - 1; i++){
        X.push_back(x1 + (x2 - x1) * (i + 1) / order);
    }
    X.push_back(x2);
}

double basis_function::evaluate(double x, const vector<double>& nodal_values){
    vector<double> shape_functions = shape_function_values(x);
    double ans = 0.0;
    for (int i = 0; i < ORDER + 1; i++)
        ans += shape_functions[i] * nodal_values[i];

    return ans;
}

vector<double> basis_function::shape_function_values(double x){
    vector<double> shape;
    shape.resize(ORDER + 1);
    for (int i = 0; i < ORDER + 1; i++)
    {
        double val = 1;
        for (int j = 0; j < ORDER + 1; j++)
        {
            if(i == j)
                continue;

            val *= (x - X[j]) / (X[i] - X[j]);
        }
        shape[i] = val;
    }

    return shape;
}

vector<double> basis_function::shape_function_gradients(double x){
    // vector<double> shape_functions = shape_function_values(x);
    vector<double> shape_grads;
    shape_grads.resize(ORDER + 1);
    for (int j = 0; j < ORDER + 1; j++)
    {
        double val = 0.0;
        for (int i = 0; i < ORDER + 1; i++){
            if(i == j)
                continue;
            double _val = 1.0 / (X[j] - X[i]);
            for (int m = 0; m < ORDER + 1; m++){
                if(m == i || m == j)
                    continue;
                _val *= (x - X[m]) / (X[j] - X[m]);
            }
            val += _val;
        }
        shape_grads[j] = val;
    }

    return shape_grads;
}

LegendreBasis::LegendreBasis(int order):ORDER(order){
    if(order < 1){
        throw logic_error("Order cannot be zero");
    }
    stringstream ss;
    ss << "octave --eval \"cd ../lib/; JacobiGL(0, 0, " << order << "); save('-ascii', 'LGLPOINTS', 'ans');\"";
    system(ss.str().c_str());
    ifstream file("../lib/LGLPOINTS");
    double x;
    LGLPoints.resize(order + 1, 0.0);
    int _i_ = 0;
    while (file >> x)
        LGLPoints[_i_++] = ((fabs(x) < 1e-14) ? 0 : x);
}

double LegendreBasis::a_val(int n, int alpha, int beta){
    return (2.0 / (2.0 * n + alpha + beta)) * sqrt((n * (n + alpha + beta) * (n + alpha) * (n + beta)) / ((2.0 * n + alpha + beta - 1.0) * (2.0 * n + alpha + beta + 1.0)));
}

double LegendreBasis::Pn(double r, int n, int alpha, int beta){
    if(n == 0)
        return sqrt(pow(2, -1 * alpha - beta - 1) * tgamma(alpha + beta + 2) / (tgamma(alpha + 1) * tgamma(beta + 1)));

    if(n == 1)
        return 0.5 * Pn(r, 0, alpha, beta) * ((alpha + beta + 2.0) * r + (alpha - beta)) * sqrt((alpha + beta + 3.0) / ((alpha + 1.0) * (beta + 1.0)));
    
    double an = a_val(n, alpha, beta), an_1 = a_val(n - 1, alpha, beta);
    double bn = (pow(beta, 2) - pow(alpha, 2)) * 1.0 / ((2 * (n - 1) + alpha + beta) * (2 * (n - 1) + alpha + beta + 2));

    return (1.0 / an) * ((r - bn) * Pn(r, n - 1, alpha, beta) - an_1 * Pn(r, n - 2, alpha, beta));
}

double LegendreBasis::GradPn(double r, int n, int alpha, int beta){
    if(n == 0)
        return 0.0;
    
    return sqrt(n * (n + alpha + beta + 1.0)) * Pn(r, n - 1, alpha + 1, beta + 1);
}

/*

inline double a_val(int n, int alpha, int beta){
    return (2.0 / (2.0 * n + alpha + beta)) * sqrt((n * (n + alpha + beta) * (n + alpha) * (n + beta)) / ((2.0 * n + alpha + beta - 1.0) * (2.0 * n + alpha + beta + 1.0)));
}

double Pn(double r, int n, int alpha, int beta){
    if(n == 0)
        return sqrt(pow(2, -1 * alpha - beta - 1) * tgamma(alpha + beta + 2) / (tgamma(alpha + 1) * tgamma(beta + 1)));

    if(n == 1)
        return 0.5 * Pn(r, 0, alpha, beta) * ((alpha + beta + 2.0) * r + (alpha - beta)) * sqrt((alpha + beta + 3.0) / ((alpha + 1.0) * (beta + 1.0)));
    
    double an = a_val(n, alpha, beta), an_1 = a_val(n - 1, alpha, beta);
    double bn = (pow(beta, 2) - pow(alpha, 2)) * 1.0 / ((2 * (n - 1) + alpha + beta) * (2 * (n - 1) + alpha + beta + 2));

    return (1.0 / an) * ((r - bn) * Pn(r, n - 1, alpha, beta) - an_1 * Pn(r, n - 2, alpha, beta));
}

double GradPn(double r, int n, int alpha, int beta){
    if(n == 0)
        return 0.0;
    
    return sqrt(n * (n + alpha + beta + 1.0)) * Pn(r, n - 1, alpha + 1, beta + 1);
}
*/