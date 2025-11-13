#include "equation.h"
#include <cmath>
#include <string>
#include <thread>
#include "DU_postprocessor/post_processor.h"

const double PI = 3.14159265358979323846;

// d^2u/dt^2 + 2 * du/dt = d^2u/dx^2 + 2 * du/dx - 3 * u
int main() {
    EquationData data;
    data.a = 1;
    data.b = 2;
    data.c = -3;
    data.d = 2;
    data.l = PI / 2;
    data.alpha = 0;
    data.beta = 1;
    data.gamma = 0;
    data.delta = 1;
    data.psi1 = [](double x) { return std::exp(-x) * std::cos(x); };
    data.psi2 = [](double x) { return -std::exp(-x) * std::cos(x); };
    data.psi1_dir1 = [](double x) { return -std::exp(-x) * std::sin(x) - std::exp(-x) * std::cos(x); };
    data.psi1_dir2 = [](double x) { return 2 * std::exp(-x) * std::sin(x); };
    data.phi0 = [](double t) { return std::exp(-t) * std::cos(2. * t); };
    data.phil = [](double t) { return 0; };
    data.f = [](double x, double t) { return 0; };
    data.exact_solution = [](double x, double t) { return std::exp(-t - x) * std::cos(x) * std::cos(2. * t); };
    data.bound_type = "a1p2";
    data.approximation = "p1";

    int meth;
    std::cin >> meth;
    std::unique_ptr<Hyperbolic> solver = std::make_unique<Hyperbolic>(data,
        (meth == 1) ? SolverType::Explicit : SolverType::Implicit);
    Hyperbolic ex_solver(data, SolverType::Exact);

    int N;
    double sigma;
    double T;
    std::cout << "Enter N, sigma, T: ";
    std::cin >> N >> sigma >> T;
    std::string method;

    std::thread solver_thread([&]() {
        method = solver->solve(N, sigma, T, 1);
    });
    std::thread exaxt_thread([&]() {
        ex_solver.solve(N, sigma, T, 1);
    })
    ;
    solver_thread.join();
    exaxt_thread.join();

    PostProcessor pp_cn("result.txt");
    pp_cn.graphics_u(data.l / N, solver->_tau, method, 6);
}