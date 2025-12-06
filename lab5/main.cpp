#include <cmath>

#include "equation.h"

int main() {
    EquationData data;

    double a = 1.;
    double c = 1.;

    // data.bound_type = Second;
    // data.exact_solution = [](double x, double t) { return std::exp(-0.5 * t) * std::sin(x); };
    // data.f = [](double x, double t) { return 0.5 * std::exp(-0.5 * t) * cos(x); };
    // data.psi = [](double x) { return std::sin(x); };
    // data.l = M_PI;
    // data.phi0 = [data](double t) { return std::exp(-0.5 * t); };
    // data.phil = [data](double t) { return -std::exp(-0.5 * t);};

    // data.bound_type = Second;
    // data.phi0 = [](double t) { return 0; };
    // data.phil = [](double t) { return 1; };
    // data.psi = [](double x) { return x + std::sin(M_PI * x); };
    // data.f = [](double x, double t) { return 0; };
    // data.l = 1.0;
    // data.exact_solution = [data, a](double x, double t) { return x + std::exp(-a * M_PI * M_PI * t) * std::sin(M_PI * x); };

    data.bound_type = Second;
    data.phi0 = [](double t) { return std::sin(t); };
    data.phil = [](double t) { return -std::sin(t); };
    data.psi = [](double x) { return 0; };
    data.f = [](double x, double t) { return std::cos(x) * (std::cos(t) + std::sin(t)); };
    data.l = M_PI / 2;
    data.exact_solution = [data, a](double x, double t) { return std::sin(t) * std::cos(x); };

    int N = 50;
    double sigma = 0.1;
    double T = 15;

    Parabolic exact_solver(data, SolverType::Exact);
    exact_solver.solve(N, sigma, a, T);

    int type;
    std::cin >> type;

    // Parabolic explicit_solver(data, SolverType::Explicit);
    // explicit_solver.solve(N, sigma, a, T);

    // ParabolicPostProcessor pp("explicit.txt");
    // pp.graphics_u(data.l / N, "explicit", 6);

    // Parabolic implicit_solver(data, SolverType::Implicit);
    // implicit_solver.solve(N, sigma, a, T);

    // ParabolicPostProcessor pp_imp("implicit.txt");
    // pp_imp.graphics_u(data.l / N, "implicit", 6);

    Parabolic cn_solver(data, SolverType(type));
    std::string method = cn_solver.solve(N, sigma, a, T);

    PostProcessor pp_cn("result.txt");
    pp_cn.graphics_u(data.l / N, cn_solver._tau, method, 6);
}