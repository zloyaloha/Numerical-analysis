#include <cstdlib>

#include "solver.h"

double f(const double& x, const double& y, const double& dy) { return std::tan(x) * dy - 2 * y; }

double exact(const double& x) { return std::sin(x) + 2 - std::sin(x) * std::log((1 + std::sin(x)) / (1 - std::sin(x))); }

int main() {
    const double a = 0;
    const double b = M_PI / 6.0;
    const double alpha = 2;
    const double beta = 2.5 - 0.5 * std::log(3);
    const double h = 0.01;

    ShootingSolver shooting_solver(a, b, alpha, beta, 0.00000001);

    std::vector<Point> res = shooting_solver.solve(f, h);
    std::vector<Point> res_h2 = shooting_solver.solve(f, h / 2);

    shooting_solver.error(res, res_h2, exact, 2);
    shooting_solver.grafics(res, exact);

    FiniteDifferenceSolver finite_difference_solver(a, b, alpha, beta, 0.001);

    std::vector<Point> res_finite_difference = finite_difference_solver.solve(f, h);
    std::vector<Point> res_finite_difference_h2 = finite_difference_solver.solve(f, h / 2);

    finite_difference_solver.error(res_finite_difference, res_finite_difference_h2, exact, 2);
    finite_difference_solver.grafics(res_finite_difference, exact);
}