#include "equation.h"
#include <cmath>

double exact(const double& x)
{
    return std::pow(x, 2) + x + 1;
}

double f(const double& x, const double& y1, const double& y2)
{
    return (2 * x * y2 - 2 * y1) / (std::pow(x, 2) - 1);
}

int main() {
    EulerSolver solver_euler(2.0, 7.0, 5.0, 3.0);
    std::vector<Point> res = solver_euler.solve(f, 0.1);
    std::vector<Point> res_h2 = solver_euler.solve(f, 0.05);

    solver_euler.error(res, res_h2, exact, 2);
    solver_euler.grafics(res, exact);

    EulerCauchySolver solver_euler_cauch(2.0, 7.0, 5.0, 3.0);
    res = solver_euler_cauch.solve(f, 0.1);
    res_h2 = solver_euler_cauch.solve(f, 0.05);

    solver_euler_cauch.error(res, res_h2, exact, 2);
    solver_euler_cauch.grafics(res, exact);

    RungeKuttaSolver solver_runge_kutta(2.0, 7.0, 5.0, 3.0);
    res = solver_runge_kutta.solve(f, 0.1);
    res_h2 = solver_runge_kutta.solve(f, 0.05);

    solver_runge_kutta.error(res, res_h2, exact, 2);
    solver_runge_kutta.grafics(res, exact);

    AdamsSolver solver_adams(2.0, 7.0, 5.0, 3.0);
    res = solver_adams.solve(f, 0.1);
    res_h2 = solver_adams.solve(f, 0.05);

    solver_adams.error(res, res_h2, exact, 2);
    solver_adams.grafics(res, exact);
}