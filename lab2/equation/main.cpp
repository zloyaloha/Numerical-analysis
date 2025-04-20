#include "equation.h"
#include <fstream>

double f(const double& x) {
    return x * std::exp(x) + std::pow(x, 2) - 1;
}

double df(const double& x) {
    return std::exp(x) + x * std::exp(x) + 2 * x;
}

double phi(const double& x) {
    return 1 / (std::exp(x) + x);
    // return std::sqrt(1 - x * std::exp(x));
    // return (1 - x * x) / std::exp(x);
}

int main() {
    double eps;
    std::cin >> eps;
    NonLinear nl(eps, 0.5);
    double res = nl.Newton(f, df);
    std::cout << res << std::endl;

    if (std::fabs(f(res)) < eps) {
        std::cout << "Проверка выполнена успешна\n";
    } else {
        std::cout << "Корень найден не верно\n";
    }

    res = nl.SimpleIter(phi);
    std::cout << res << std::endl;

    if (fabs(f(res)) < eps) {
        std::cout << "Проверка выполнена успешна\n";
    } else {
        std::cout << "Корень найден не верно " << f(res) << std::endl;
    }
}