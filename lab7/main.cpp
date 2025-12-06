#include <cmath>
#include <memory>
#include <string>

#include "DU_postprocessor/post_processor.h"
#include "equation.h"

int main() {
    EquationData eqData;
    eqData.phi0 = [](double y) { return std::cos(y); };
    eqData.phi1 = [](double y) { return 0; };
    eqData.phi2 = [](double x) { return std::cos(x); };
    eqData.phi3 = [](double x) { return 0; };
    eqData.exact = [](double x, double y) { return std::cos(x) * std::cos(y); };
    eqData.l = M_PI / 2;

    int meth;
    std::cout << "Enter method (1 - Seidel, 2 - Leibmann, 3 - SOR): ";
    std::cin >> meth;
    std::unique_ptr<Elliptic> solver = std::make_unique<Elliptic>(eqData, SolverType(meth));
    Elliptic ex_solver(eqData, SolverType::Exact);

    int N;
    double epsilon;
    std::cout << "Enter N, eps: ";
    std::cin >> N >> epsilon;
    std::string method;

    method = solver->solve(N, epsilon);
    ex_solver.solve(N, epsilon);

    PostProcessor pp_cn("result.txt");
    pp_cn.graphics_elliptic(eqData.l / N, method);
}