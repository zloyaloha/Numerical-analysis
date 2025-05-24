#include "system.h"
#include <cmath>
#include <fstream>

double f1(const double& x1, const double& x2)
{
    return 2 * x1 - std::cos(x2);
}

double f2(const double& x1, const double& x2)
{
    return 2 * x2 - std::exp(x1);
}

double df1_dx1(const double& x1, const double& x2)
{
    return 2;
}

double df1_dx2(const double& x1, const double& x2)
{
    return sin(x2);
}

double df2_dx1(const double& x1, const double& x2)
{
    return -std::exp(x1);
}

double df2_dx2(const double& x1, const double& x2)
{
    return 2;
}

double phi1(const double& x1, const double& x2) {
    return std::cos(x2) / 2;
}

double phi2(const double& x1, const double& x2) {
    return std::exp(x1) / 2;
}

double dphi1_dx1(const double& x1, const double& x2) {
    return 0;
}

double dphi1_dx2(const double& x1, const double& x2) {
    return -std::sin(x2) / 2;
}

double dphi2_dx1(const double& x1, const double& x2) {
    return std::exp(x1) / 2;
}

double dphi2_dx2(const double& x1, const double& x2) {
    return 0;
}

int main() {
    double eps;
    std::cin >> eps;
    Matrix<std::function<double(const double&, const double&)>> J(2, 2);
    J.at(0, 0) = df1_dx1;
    J.at(0, 1) = df1_dx2;
    J.at(1, 0) = df2_dx1;
    J.at(1, 1) = df2_dx2;

    NonLinearSystem system(eps, {-0.5, -0.5});
    std::vector<double> res = system.Newton(J, f1, f2);
    if (f1(res[0], res[1]) > eps || f2(res[0], res[1]) > eps) {
        std::cerr << "Решение не прошло проверку" << std::endl;
    } else {
        std::cerr << "Решение прошло проверку" << std::endl;
    }
    for (const double& x: res) {
        std::cout << x << ' ';
    }
    std::cout << std::endl;
    
    J.at(0, 0) = dphi1_dx1;
    J.at(0, 1) = dphi1_dx2;
    J.at(1, 0) = dphi2_dx1;
    J.at(1, 1) = dphi2_dx2;

    res = system.SimpleIter(phi1, phi2, J, 0, 1, 0, 1);

    if (f1(res[0], res[1]) > eps || f2(res[0], res[1]) > eps) {
        std::cerr << "Решение не прошло проверку" << std::endl;
    } else {
        std::cerr << "Решение прошло проверку" << std::endl;
    }
    for (const double& x: res) {
        std::cout << x << ' ';
    }
    std::cout << std::endl;

}