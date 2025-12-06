#include <cmath>
#include <iostream>

#include "equation.h"

int main() {
    EquationData data;
    data.l1 = 1.0;
    data.l2 = 1.0;

    data.exact = [](double x, double y, double t) { return x * y * std::cos(t); };

    data.f = [](double x, double y, double t) { return -x * y * std::sin(t); };
    data.psi = [](double x, double y) { return x * y; };

    data.phi0 = [](double y, double t) { return 0.0; };              // x=0
    data.phi1 = [](double y, double t) { return y * std::cos(t); };  // x=1
    data.phi2 = [](double x, double t) { return 0.0; };              // y=0
    data.phi3 = [](double x, double t) { return x * std::cos(t); };  // y=1

    int N = 40;
    double T_end = 2 * M_PI;

    try {
        Parabolic2D solver_exp(data, SolverType::FractSteps, 15);
        solver_exp.solve(N, N, 0.25, T_end);
    } catch (const std::exception& e) {
        std::cerr << "Explicit Solver Error: " << e.what() << std::endl;
    }

    try {
        Parabolic2D solver_imp(data, SolverType::AlterDirection, 15);
        solver_imp.solve(N, N, 0.5, T_end);
    } catch (const std::exception& e) {
        std::cerr << "Implicit Solver Error: " << e.what() << std::endl;
    }

    try {
        Parabolic2D solver_imp(data, SolverType::Exact);
        solver_imp.solve(N, N, 0.5, T_end);
    } catch (const std::exception& e) {
        std::cerr << "Exact Solver Error: " << e.what() << std::endl;
    }

    return 0;
}