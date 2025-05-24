#include "equation.h"
#include <fstream>

double f(const double& x) {
    return x * std::exp(x) + std::pow(x, 2) - 1;
}

double df(const double& x) {1 / (std::exp(x) + x);
    return std::exp(x) + x * std::exp(x) + 2 * x;
}

double ddf(const double& x) {
    return 2 * std::exp(x) + x * std::exp(x) + 2;
}

double phi(const double& x) {
    return 1 / (std::exp(x) + x);
}

double dphi(const double& x) {
    return -(std::exp(x) + 1) / std::pow(std::exp(x) + x, 2);
}

int main() {
    double eps;
    std::cin >> eps;
    NonLinear nl(eps);
    double res = nl.Newton(f, df, ddf, 0, 0.5);
    std::cout << res << std::endl;

    if (std::fabs(f(res)) < eps) {
        std::cout << "Проверка выполнена успешна\n";
    } else {
        std::cout << "Корень найден не верно\n";
    }

    res = nl.SimpleIter(phi, dphi, 0.3, 0.5);
    std::cout << res << std::endl;

    if (fabs(f(res)) < eps) {
        std::cout << "Проверка выполнена успешна\n";
    } else {
        std::cout << "Корень найден не верно " << f(res) << std::endl;
    }
}