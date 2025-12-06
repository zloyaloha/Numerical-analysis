#include "equation.h"
#include "matrix.h"
#include "tridiagonal-matrix.h"

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

class FiniteDifferenceSolver : public Solver {
public:
    FiniteDifferenceSolver(const double& a, const double& b, const double& alpha, const double& beta, const double& eps);
    std::vector<Point> solve(std::function<double(const double&, const double&, const double&)> f, const double& h) override;

private:
    double a, b;
    double alpha, beta;
    double h;
    double eps;
};

ShootingSolver::ShootingSolver(const double& a_, const double& b_, const double& alpha_, const double& beta_, const double& eps_)
    : a(a_), b(b_), alpha(alpha_), beta(beta_), eps(eps_) {
    name = "Shooting";
    p = 2.;
}

FiniteDifferenceSolver::FiniteDifferenceSolver(const double& a_, const double& b_, const double& alpha_, const double& beta_,
                                               const double& eps_)
    : a(a_), b(b_), alpha(alpha_), beta(beta_), eps(eps_) {
    name = "Finite Difference";
    p = 2.;
}

std::vector<Point> ShootingSolver::solve(std::function<double(const double&, const double&, const double&)> f, const double& h) {
    auto integrate = [&](double s) -> std::vector<Point> {
        RungeKuttaSolver solver(a, alpha, s, b);
        return solver.solve(f, h);
    };

    double n_prev = 1;
    double n = -1;

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

double ShootingSolver::phi(const std::vector<Point>& sol) { return sol.back().y1 - beta; }

double ShootingSolver::secant(double s_prev, double s, const std::vector<Point>& sol_prev, const std::vector<Point>& sol) {
    double phi_prev = phi(sol_prev);
    double phi_curr = phi(sol);
    return s - (s - s_prev) / (phi_curr - phi_prev) * phi_curr;
}

std::vector<Point> FiniteDifferenceSolver::solve(std::function<double(const double&, const double&, const double&)> f,
                                                 const double& h) {
    int N = (b - a) / h;
    Tridiagonal tr;

    std::vector<double> d(N - 1, 0.0);     // правая часть: r(x) = 0
    std::vector<double> lower(N - 2);      // поддиагональ
    std::vector<double> main_diag(N - 1);  // главная диагональ
    std::vector<double> upper(N - 2);      // наддиагональ

    for (int i = 1; i < N; ++i) {
        double x_i = a + i * h;

        double p = -std::tan(x_i);
        double q = 2.0;
        double r = 0.0;

        double A = 1.0 / (h * h) - p / (2.0 * h);
        double B = -2.0 / (h * h) + q;
        double C = 1.0 / (h * h) + p / (2.0 * h);

        if (i > 1) lower[i - 2] = A;
        main_diag[i - 1] = B;
        if (i < N - 1) upper[i - 1] = C;

        d[i - 1] = r;

        if (i == 1) d[i - 1] -= A * alpha;
        if (i == N - 1) d[i - 1] -= C * beta;
    }

    Matrix<double> A(main_diag, lower, upper, d);
    std::vector<double> y_inner = tr.solve(A);

    std::vector<Point> solution;
    solution.push_back({a, alpha, 0.0});
    for (int i = 1; i < N; ++i) {
        double x_i = a + i * h;
        solution.push_back({x_i, y_inner[i - 1], 0.0});
    }
    solution.push_back({b, beta, 0.0});

    return solution;
}