#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <functional>
#include <stdexcept>

namespace {
    const int MAX_ITER = 100;
}

class NonLinear {
private:
    double _eps;
    double _x0;

public:
    NonLinear(const double& eps, const double& x0): _eps(eps), _x0(x0) {}

    double Newton(std::function<double(const double&)> f, std::function<double(const double&)> df)
    {
        double x = _x0;
        for (int i = 0; i < MAX_ITER; ++i) {
            double fx = f(x);
            double dfx = df(x);

            double delta = fx / dfx;
            x -= delta;

            if (fabs(delta) < _eps) {
                std::cout << "Решение найдено за " << i << " итераций\n";
                return x;
            }
        }
        std::cerr << "Решение не найдено за " << MAX_ITER << " итераций\n";
        return x;
    }

    // на интервале от (0.4; 0.6) модуль производной меньше 0.7
    double SimpleIter(std::function<double(const double&)> phi)
    {
        double x = _x0;
        double x_next = 0;
        for (int i = 0; i < MAX_ITER; ++i) {
            x_next = phi(x);
            if ((0.7 / (1 - 0.7)) * fabs(x_next - x) < _eps) {
                std::cout << "Решение найдено за " << i + 1 << " итераций." << std::endl;
                return x_next;
            }
            
            x = x_next;
        }
        std::cerr << "Достигнуто максимальное число итераций!" << std::endl;
        return x;
    }
};