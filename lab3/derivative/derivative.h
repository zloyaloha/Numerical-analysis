#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

class Derivative {
public:
    double findFirstDerivative(const std::vector<double>& x, const std::vector<double>& y, double x_star) {
        int idx = -1;
        for (int i = 0; i < x.size() - 1; ++i) {
            if (x[i] < x_star && x_star <= x[i + 1]) {
                idx = i + 1;
                break;
            }
        }

        if (x[0] == x_star) {
            idx = 0;
        }

        if (x[idx] == x_star) {
            double right, left;
            if (idx + 1 < x.size()) {
                right = (y[idx + 1] - y[idx]) / (x[idx + 1] - x[idx]);
                std::cout << "Правосторонняя первая производная " << right << '\n';
            }
            if (idx != 0) {
                left = (y[idx] - y[idx - 1]) / (x[idx] - x[idx - 1]);
                std::cout << "Левосторонняя первая производная " << left << '\n';
            }
        }

        if (idx == -1) {
            throw std::out_of_range("x* вне диапазона x.");
        }

        if (idx > 0 && idx < x.size() - 2) {
            double drob1 = (y[idx + 1] - y[idx]) / (x[idx + 1] - x[idx]);
            double drob21 = (y[idx + 2] - y[idx + 1]) / (x[idx + 2] - x[idx + 1]);
            double drob22 = (y[idx + 1] - y[idx]) / (x[idx + 1] - x[idx]);
            double drob3 = (2 * x_star - x[idx] - x[idx + 1]) / (x[idx + 2] - x[idx]);

            return drob1 + (drob21 - drob22) * drob3;
        } else {
            throw std::out_of_range("x* слишком с краю");
        }
    }

    double findSecondDerivative(const std::vector<double>& x, const std::vector<double>& y, double x_star) {
        int idx = -1;
        for (int i = 0; i < x.size() - 1; ++i) {
            if (x[i] <= x_star && x_star <= x[i + 1]) {
                idx = i;
                break;
            }
        }

        if (idx > 0 && idx < x.size() - 2) {
            double drob11 = (y[idx + 2] - y[idx + 1]) / (x[idx + 2] - x[idx + 1]);
            double drob12 = (y[idx + 1] - y[idx]) / (x[idx + 1] - x[idx]);
            double drob2 = 2 / (x[idx + 2] - x[idx]);
            return (drob11 - drob12) * drob2;
        } else {
            throw std::out_of_range("x* слишком с краю");
        }
    }
};