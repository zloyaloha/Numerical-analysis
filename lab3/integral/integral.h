#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>
#include <cmath>

class Integral {
private:
    double a, b;
public:
    Integral(double a_, double b_) : a(a_), b(b_) {}

    double findIntegralRectangle(std::function<double(double)> f, double h)
    {
        double sum = 0;
        for (double x = a; x < b; x += h) {
            double xi = (x + x + h) / 2;
            sum += f(xi);
        }
        return sum * h;
    }

    double findIntegralTrapezoid(std::function<double(double)> f, double h)
    {
        double sum = (f(a) + f(b)) / 2.0;
        for (double x = a + h; x < b; x += h) {
            sum += f(x);
        }
        return sum * h;
    }

    double findIntegralSimpson(std::function<double(double)> f, double h)
    {
        int n = (b - a) / h;
        if (n % 2 != 0) ++n; // делаем n чётным
        h = (b - a) / n;

        double sum = f(a) + f(b);
        for (int i = 1; i < n; ++i) {
            double x = a + i * h;
            sum += f(x) * (i % 2 == 0 ? 2 : 4);
        }
        return sum * h / 3.0;
    }

    double rungeRomberg(double I_h, double I_h2, int p) {
        return I_h2 + (I_h2 - I_h) / (std::pow(2, p) - 1);
    }
};  