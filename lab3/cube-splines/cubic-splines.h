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
public:
    CubicSplines(const std::vector<double>& x, const std::vector<double>& y) {
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
};