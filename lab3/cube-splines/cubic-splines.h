#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>
#include "tridiagonal-matrix.h"

class CubicSplines {
private:
    struct SplineSegment {
        double a, b, c, d;
        double x1, x2;
    };
    std::vector<SplineSegment> segments;

    std::vector<double> x_points, y_points;
    double eps{0.01};

    void validateInput(const std::vector<double>& x, const std::vector<double>& y) {
        if (x.size() != y.size()) {
            throw std::invalid_argument("x and y vectors must have the same size");
        }
        if (x.size() < 2) {
            throw std::invalid_argument("At least 2 points are required for spline interpolation");
        }

        for (size_t i = 1; i < x.size(); ++i) {
            if (x[i] <= x[i-1]) {
                throw std::invalid_argument("x values must be strictly increasing");
            }
        }
    }

    bool checkEdgeNodes() const {
        if (std::abs(y_points[0] - evaluate(x_points[0])) >= eps) {
            return false;
        }

        if (std::abs(y_points.back() - evaluate(x_points.back())) >= eps) {
            return false;
        }

        return true;
    }

    bool checkInternalNodesContinuity() const {
        for (size_t i = 1; i < x_points.size() - 1; ++i) {
            if (std::abs(y_points[i] - evaluate(x_points[i])) >= eps) {
                return false;
            }
        }
        return true;
    }

    bool checkFirstDerivativeContinuity() const {
        for (size_t i = 1; i < x_points.size() - 1; ++i) {
            double left_deriv = derivative(x_points[i], i-1, 1);
            double right_deriv = derivative(x_points[i], i, 1);

            if (std::abs(left_deriv - right_deriv) >= eps) {
                return false;
            }
        }
        return true;
    }

    bool checkSecondDerivativeContinuity() const {
        for (size_t i = 1; i < x_points.size() - 1; ++i) {
            double left_deriv = derivative(x_points[i], i-1, 2);
            double right_deriv = derivative(x_points[i], i, 2);

            if (std::abs(left_deriv - right_deriv) >= eps) {
                return false;
            }
        }
        return true;
    }

    bool checkEdgeSecondDerivatives() const {
        if (std::abs(derivative(x_points[0], 0, 2)) >= eps) {
            return false;
        }

        if (std::abs(derivative(x_points.back(), segments.size()-1, 2)) >= eps) {
            return false;
        }

        return true;
    }

    double derivative(double x, size_t segment_idx, int order) const {
        if (segment_idx >= segments.size()) {
            throw std::out_of_range("Invalid segment index");
        }

        const SplineSegment& s = segments[segment_idx];
        double dx = x - s.x1;

        switch (order) {
            case 1:  // First derivative
                return s.b + 2*s.c*dx + 3*s.d*dx*dx;
            case 2:  // Second derivative
                return 2*s.c + 6*s.d*dx;
            default:
                throw std::invalid_argument("Only 1st and 2nd derivatives are supported");
        }
    }

public:

    CubicSplines(const std::vector<double>& x, const std::vector<double>& y) : x_points(x), y_points(y) {
        validateInput(x, y);
        const int n = x.size() - 1;
        std::vector<double> h(n);
        for (int i = 0; i < n; ++i) {
            h[i] = x[i + 1] - x[i];
        }

        Matrix<double> system(n - 1, n, 0.0); // (n-1)x(n) — расширенная матрица
        for (int i = 0; i < n - 1; ++i) {
            system.at(i, i)     = 2.0 * (h[i] + h[i + 1]);      // диагональ
            if (i > 0) system.at(i, i - 1) = h[i];              // нижняя диагональ
            if (i < n - 2) system.at(i, i + 1) = h[i + 1];      // верхняя диагональ

            system.at(i, n - 1) = 3.0 * ((y[i + 2] - y[i + 1]) / h[i + 1] -
                                        (y[i + 1] - y[i]) / h[i]);
        }

        Tridiagonal tr;
        std::vector<double> internal_c = tr.solve(system);

        std::vector<double> c(n + 1, 0.0);
        for (int i = 1; i < n; ++i) {
            c[i] = internal_c[i - 1];
        }

        segments.resize(n);
        for (int i = 0; i < n; ++i) {
            segments[i].x1 = x[i];
            segments[i].x2 = x[i + 1];
            segments[i].a = y[i];
            segments[i].c = c[i];
            segments[i].d = (c[i + 1] - c[i]) / (3.0 * h[i]);
            segments[i].b = (y[i + 1] - y[i]) / h[i] - h[i] * (2.0 * c[i] + c[i + 1]) / 3.0;
        }
    }

    double evaluate(double x_star) const {
        // std::cout << segments[segments.size() - 1].x1 << ' ' << x_star << ' ' << segments[0].x1 << std::endl;
        // if (x_star > segments[segments.size() - 1].x1 || x_star < segments[0].x1) {
        //     throw std::range_error("X not in range");
        // }
        auto it = std::upper_bound(segments.begin(), segments.end(), x_star,
            [](double val, const SplineSegment& seg) { return val < seg.x1; });

        if (it == segments.begin()) return segments.front().a;

        if (it == segments.end()) {
            const auto& last = segments.back();
            double dx = x_star - last.x1;
            return last.a + last.b*dx + last.c*dx*dx + last.d*dx*dx*dx;
        }

        const SplineSegment& s = *(it - 1);
        double dx = x_star - s.x1;
        return s.a + s.b*dx + s.c*dx*dx + s.d*dx*dx*dx;
    }

    void printPolynomials() const {
        std::cout << "Кусочно-кубические многочлены сплайна:\n";
        std::cout << "Интервал [x_i, x_{i+1}]:\n";
        std::cout << "S_i(x) = a + b*(x - x_i) + c*(x - x_i)^2 + d*(x - x_i)^3\n\n";

        for (size_t i = 0; i < segments.size(); ++i) {
            const auto& s = segments[i];
            std::cout << "Интервал [" << s.x1 << ", " << s.x2 << "]:\n";
            std::cout << "  a = " << s.a << "\n";
            std::cout << "  b = " << s.b << "\n";
            std::cout << "  c = " << s.c << "\n";
            std::cout << "  d = " << s.d << "\n";
            std::cout << "  Многочлен: S_" << i << "(x) = "
                    << s.a << " + "
                    << s.b << "*(x - " << s.x1 << ") + "
                    << s.c << "*(x - " << s.x1 << ")^2 + "
                    << s.d << "*(x - " << s.x1 << ")^3\n\n";
        }
    }

    void writeToFile(const std::string& filename, double step = 0.01) const {
        std::ofstream file(filename);
        for (double xi = segments.front().x1; xi <= segments.back().x2; xi += step) {
            file << xi << " " << evaluate(xi) << "\n";
        }
    }

    bool validateSpline() const {
        return checkEdgeNodes() &&
               checkInternalNodesContinuity() &&
               checkFirstDerivativeContinuity() &&
               checkSecondDerivativeContinuity() &&
               checkEdgeSecondDerivatives();
    }
};