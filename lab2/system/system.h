#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <functional>
#include <stdexcept>
#include "matrix.h"
#include "lu.h"

namespace {
    const int MAX_ITER = 100;
}

class NonLinearSystem {
private:
    double _eps;
    const std::vector<double> _x0;

public:
    NonLinearSystem(const double& eps, const std::vector<double> x0): _eps(eps), _x0(x0) {}

    std::vector<double>  Newton(Matrix<std::function<double(const double&, const double&)>> J, std::function<double(const double&, const double&)> f1, std::function<double(const double&, const double&)> f2)
    {
        LU lu;
        double x1 = _x0[0], x2 = _x0[1];
        Matrix<double> system(2, 3);
        int iter = 0;
        while (iter < MAX_ITER) {
            system.at(0, 0) = J.at(0, 0)(x1, x2);
            system.at(0, 1) = J.at(0, 1)(x1, x2);
            system.at(0, 2) = -f1(x1, x2);

            system.at(1, 0) = J.at(1, 0)(x1, x2);
            system.at(1, 1) = J.at(1, 1)(x1, x2);
            system.at(1, 2) = -f2(x1 , x2);

            double norm = std::sqrt(std::pow(system.at(0, 2), 2) + std::pow(system.at(1, 2), 2));

            if (norm < _eps) {
                std::cout << "Решение найдено за " << iter << " итераций" << std::endl;
                return {x1, x2};
            }

            std::vector<double> x = lu.solve(system);

            x1 += x[0]; x2 += x[1];
            ++iter;
        }
        std::cout << "Решение не удалось найти за " << MAX_ITER << " итераций" << std::endl;
        return {x1, x2};
    }

    std::vector<double> SimpleIter(
        std::function<double(const double&, const double&)> phi1, 
        std::function<double(const double&, const double&)> phi2, 
        Matrix<std::function<double(const double&, const double&)>> J, 
        double A1, double A2, double B1, double B2)
        {
        double x1 = _x0[0], x2 = _x0[1];
        double max_q = 0;
        for (double i = A1; i <= A1; i += 0.05) {
            for (double j = A2; j <= A2; j += 0.05) {
                for (int row = 0; row < J.getRows(); ++row) {
                    double sum = 0;
                    for (int col = 0; col < J.getCols(); ++col) {
                        sum += std::fabs(J.at(row, col)(i, j));
                    }
                    if (max_q < sum) {
                        max_q = sum;
                        x1 = i;
                        x2 = j;
                    }
                }
            }
        }
        std::cout << "q = " << max_q << std::endl;
        if (max_q >= 1) {
            std::cerr << "Может не сойтись" << std::endl;
            std::cout << max_q << std::endl;
        }

        int iter = 0;
        while (iter < MAX_ITER) {
            double x1_new = phi1(x1, x2);
            double x2_new = phi2(x1, x2);

            double delta = std::max(std::fabs(x1_new - x1), std::fabs(x2_new - x2));

            if (delta < _eps) {
                std::cout << "Решение найдено за " << iter << " итераций" << std::endl;
                return {x1, x2};
            }

            x1 = x1_new;
            x2 = x2_new;
            iter++;
        }
        std::cout << "Решение не удалось найти за " << MAX_ITER << " итераций" << std::endl;
        return {x1, x2};
    }
};