#include <iostream>
#include <iterator>
#include <vector>
#include <cmath>
#include <algorithm>
#include "newton_lagrange.h"

double f(const double& x)
{
    return std::log(x) + x;
}

int main()
{
    std::vector<double> X1 = {0.1, 0.5, 0.9, 1.3};
    std::vector<double> X2 = {0.1, 0.5, 1.1, 1.3};
    LagrangeInterpolator inter(f, X1);

    std::vector<double> Y1(X1.size());
    std::vector<double> Y2(X2.size());

    std::transform(X1.begin(), X1.end(), Y1.begin(), f);
    std::transform(X2.begin(), X2.end(), Y2.begin(), f);

    double x = 0.8;

    double val = inter.interpolate(x);
    double real_val = f(x);

    std::cout << "Приближение в точке x = " << x << " методом Лагранжа P(x) = " << val << '\n';
    std::cout  << "Реальное значение функции в точке f(x) = " << real_val << std::endl;
    std::cout << "Погрешность " << std::fabs(val - real_val) << std::endl;

    inter.setNewPoints(X2);
    val = inter.interpolate(x);

    std::cout << "Приближение в точке x = " << x << " методом Лагранжа P(x) = " << val << '\n';
    std::cout  << "Реальное значение функции в точке f(x) = " << real_val << std::endl;
    std::cout << "Погрешность " << std::fabs(val - real_val) << std::endl;

    // ==================

    NewtonInterpolator newton_inter(f, X1);

    val = newton_inter.interpolate(x);

    std::cout << "Приближение в точке x = " << x << " методом Ньютона P(x) = " << val << '\n';
    std::cout  << "Реальное значение функции в точке f(x) = " << real_val << std::endl;
    std::cout << "Погрешность " << std::fabs(val - real_val) << std::endl;

    newton_inter.setNewPoints(X2);
    val = newton_inter.interpolate(x);

    std::cout << "Приближение в точке x = " << x << " методом Ньютона P(x) = " << val << '\n';
    std::cout  << "Реальное значение функции в точке f(x) = " << real_val << std::endl;
    std::cout << "Погрешность " << std::fabs(val - real_val) << std::endl;
}