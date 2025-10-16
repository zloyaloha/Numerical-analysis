#include <iostream>
#include "matrix.h"
#include <vector>
#include <stdexcept>

class Tridiagonal {
public:
    Tridiagonal() = default;

    template <class T>
    typename std::enable_if<std::is_arithmetic<T>::value, std::vector<double>>::type
    solve(Matrix<T>& A)
    {
        int n = A.getRows();
        double a = 0, b = A.at(0, 0), c = A.at(0, 1), d = A.at(0, A.getCols() - 1);
        std::vector<double> P(n), Q(n);
        P[n - 1] = 0;
        P[0] = -c / b; Q[0] = d / b;
        int counter = 0;
        for (int i = 1; i < n; ++i) {
            a = A.at(i, i - 1);
            b = A.at(i, i);
            c = A.at(i, i + 1);
            d = A.at(i, A.getCols() - 1);

            if (a == 0 || c == 0) {
                std::cerr << "the method may not converge. a or c is zero" << std::endl;
            }
            if (b < a + c) {
                std::cerr << "the method may not converge. b is lower then a + c" << std::endl;
            }
            if (b > a + c) {
                ++counter;
            }

            P[i] = -c / (b + a * P[i - 1]);
            Q[i] = (d - a * Q[i - 1]) / (b + a * P[i - 1]);
        }
        if (counter == 0) {
            std::cerr << "the method may not converge. b is not strongly bigger then a + c at least once" << std::endl;
        }
        std::vector<double> x(n);
        x[n - 1] = Q[n - 1];
        for (int i = n - 2; i >= 0; --i) {
            x[i] = P[i] * x[i + 1] + Q[i];
        }
        return x;
    }
};

inline std::vector<double> tma(const std::vector<double>& a,
                   const std::vector<double>& b,
                   const std::vector<double>& c,
                   const std::vector<double>& d)
{
    int n = b.size();
    std::vector<double> alpha(n), beta(n), x(n);

    alpha[0] = -c[0] / b[0];
    beta[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        double denom = b[i] + a[i] * alpha[i - 1];
        alpha[i] = -c[i] / denom;
        beta[i] = (d[i] - a[i] * beta[i - 1]) / denom;
    }

    x[n - 1] = beta[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }

    return x;
}
