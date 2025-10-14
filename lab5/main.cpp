#include "equation.h"
#include <cmath>

int main() {
    EquationData data;

    data.bound_type = Second;
    data.exact_solution = [](double x, double t) { return std::exp(-0.5 * t) * std::sin(x); };
    data.f = [](double x, double t) { return 0.5 * std::exp(-0.5 * t) * cos(x); };
    data.phi0 = [](double t) { return -exp(-0.5 * t); };
    data.phil = [](double t) { return -exp(-0.5 * t); };
    data.psi = [](double x) { return std::sin(x); };
    data.l = M_PI;

    Parabolic explicit_solver(data, SolverType::Explicit);
    auto res_ex = explicit_solver.solve(20, 100, 1);
    explicit_solver.grafics_u(res_ex, 1, "explicit");

    Parabolic implicit_solver(data, SolverType::Implicit);
    auto res_im = implicit_solver.solve(20, 100, 1);
    implicit_solver.grafics_u(res_im, 1, "implicit");

    std::cout << "u(x=0.5, t=1.0) ≈ " << res_ex.back()[5] << std::endl;
    std::cout << "u(x=0.5, t=1.0) ≈ " << res_im.back()[5] << std::endl;
    std::cout << "analytic: " << data.exact_solution(0.5, 1.0) << std::endl;

    Parabolic exact_solver(data, SolverType::Exact);
    auto res_exact = exact_solver.solve(20, 100, 1);
    exact_solver.grafics_u(res_exact, 1, "exact");
}