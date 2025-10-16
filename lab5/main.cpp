#include "equation.h"
#include <cmath>

int main() {
    EquationData data;

    double a = 1.;
    double c = 1.;

    data.bound_type = Second;
    data.exact_solution = [](double x, double t) { return std::exp(-0.5 * t) * std::sin(x); };
    data.f = [](double x, double t) { return 0.5 * std::exp(-0.5 * t) * cos(x); };
    data.phi0 = [data](double t) { return data.exact_solution(t, 0); };
    data.psi = [](double x) { return std::sin(x); };
    data.l = M_PI;
    data.phil = [data](double t) { return data.exact_solution(t, data.l); };

    int N = 100;
    double sigma = 0.4;
    double T = 0.5;

    Parabolic exact_solver(data, SolverType::Exact);
    exact_solver.solve(N, sigma, a, T);

    Parabolic explicit_solver(data, SolverType::Explicit);
    explicit_solver.solve(N, sigma, a, T);

    ParabolicPostProcessor pp("explicit.txt");
    pp.graphics_u(data.l / N, "explicit", 6);

    Parabolic implicit_solver(data, SolverType::Implicit);
    implicit_solver.solve(N, sigma, a, T);

    ParabolicPostProcessor pp_imp("implicit.txt");
    pp_imp.graphics_u(data.l / N, "implicit", 6);

    Parabolic cn_solver(data, SolverType::Crank_Nicholson);
    cn_solver.solve(N, sigma, a, T);

    ParabolicPostProcessor pp_cn("crank-nicolson.txt");
    pp_cn.graphics_u(data.l / N, "crank_nicolson", 6);
}