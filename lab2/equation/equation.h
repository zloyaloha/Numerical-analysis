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
    NonLinear(const double& eps): _eps(eps) {}

    double Newton(std::function<double(const double&)> f, std::function<double(const double&)> df, std::function<double(const double&)> ddf, const double& a, const double& b)
    {
        if (f(a) * f(b) > 0) {
            std::cerr << "Корня на отрезке нет" << std::endl;
            return 0;
        }

        if (df(_x0) == 0) {
            std::cerr << "Производная ноль" << std::endl;
            return 0;
        }

        bool flag = false;
        for (double i = a; i < b; i += 0.05) {
            if (f(i) * ddf(i) > 0) {
                _x0 = i;
                flag = true;
            }
        }
        if (!flag) {
            std::cerr << "Может не сойтись" << std::endl;
            _x0 = a;
        }

        double x = _x0;
        for (int i = 0; i < MAX_ITER; ++i) {
            double fx = f(x);
            double dfx = df(x);

            if (dfx == 0) {
                throw std::runtime_error("diff is zero");
            }
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

    double SimpleIter(std::function<double(const double&)> phi, std::function<double(const double&)> dphi, const double& a, const double& b)
    {
        bool flag = false;
        double q = 0;
        double _x0 = 0;
        for (double i = a; i <= b; i += 0.05) {
            std::cout << std::fabs(dphi(i)) << ' ' << i << std::endl;
            if (q < std::fabs(dphi(i))) {
                q = std::fabs(dphi(i));
                _x0 = i;
            }
        }

        std::cout << "q = " << q << ' ' << "x0 = " << _x0 << std::endl;
        if (q >= 1) {
            std::cerr << "Может не сойтись" << std::endl;
        }

        double x = _x0;
        double x_next;
        for (int i = 0; i < MAX_ITER; ++i) {
            x_next = phi(x);
            if ((q / (1 - q)) * fabs(x_next - x) < _eps) {
                std::cout << "Решение найдено за " << i + 1 << " итераций." << std::endl;
                return x_next;
            }
            
            x = x_next;
        }
        std::cerr << "Достигнуто максимальное число итераций!" << std::endl;
        return x;
    }
};