#include "equation.h"

class ShootingSolver : public Solver {
public:
    ShootingSolver(const double& a, const double& b, const double& alpha, const double& beta, const double& eps);
    std::vector<Point> solve(std::function<double(const double&, const double&, const double&)> f, const double& h) override;
private:
    double a, b;
    double alpha, beta;
    double h;
    double eps;
private:
    double phi(const std::vector<Point>& sol);
    double secant(double s_prev, double s, const std::vector<Point>& sol_prev, const std::vector<Point>& sol);
};

ShootingSolver::ShootingSolver(const double& a_, const double& b_, const double& alpha_, const double& beta_, const double& eps_) : a(a_), b(b_), alpha(alpha_), beta(beta_), eps(eps_)
{
    name = "Shooting";
    p = 4.;
}

std::vector<Point> ShootingSolver::solve(std::function<double(const double&, const double&, const double&)> f, const double& h)
{
    auto integrate = [&](double s) -> std::vector<Point> {
        RungeKuttaSolver solver(a, alpha, s, b);
        return solver.solve(f, h);
    };

    double n_prev = 0.5;
    double n = -0.5;

    auto sol_prev = integrate(n_prev);
    auto sol = integrate(n);

    while (std::abs(phi(sol)) > eps) {
        double n_next = secant(n_prev, n, sol_prev, sol);
        n_prev = n;
        sol_prev = sol;
        n = n_next;
        sol = integrate(n);
    }

    return sol;
}

double ShootingSolver::phi(const std::vector<Point>& sol) {
    return sol.back().y1 - beta;
}

double ShootingSolver::secant(double s_prev, double s, const std::vector<Point>& sol_prev, const std::vector<Point>& sol) {
    double phi_prev = phi(sol_prev);
    double phi_curr = phi(sol);
    return s - (s - s_prev) / (phi_curr - phi_prev) * phi_curr;
}