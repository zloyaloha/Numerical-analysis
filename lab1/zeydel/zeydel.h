#include <iostream>
#include "matrix.h"
#include <vector>
#include <stdexcept>
#include <cmath>

const double EPS = 1e-6; // Точность вычислений
const double max_iter = 10; // Точность вычислений

class Zeydel {
public:
    Zeydel() = default;

    template <class T>
    typename std::enable_if<std::is_arithmetic<T>::value, std::vector<double>>::type
    solve(Matrix<T>& A)
    {
        int k = 0;
        double max_norm = 0;
        std::vector<double> x(A.getRows(), 0), x_prev(A.getRows(), 0);
        Matrix<double> A_(A.getRows(), A.getCols(), 0.);
        for (int i = 0; i < A.getRows(); ++i) {
            while (A.at(i, i) == 0) {
                ++k;
                A.swapRows(i, i + k);
            }
            A_.at(i, 0) = A.at(i, A.getCols() - 1) / A.at(i, i); // свободный член на последнем месте
            k = 0;
            for (int j = 0; j < A.getCols() - 1; ++j) {
                A_.at(i, j + 1) = A.at(i, j) / A.at(i, i); // тут все как обычно
            }
            A_.at(i, i + 1) = 0;
            x[i] = A_.at(i, 0);
        }

        for (int iter = 0; iter < max_iter; ++iter) {
            x_prev = x; max_norm = 0;
            for (int i = 0; i < A.getRows(); ++i) {
                double sum = A_.at(i, 0);
                for (int j = 0; j < A.getCols() - 1; ++j) {
                    sum -= A_.at(i, j + 1) * x[j];
                }
                x[i] = sum;
                max_norm = std::max(max_norm, fabs(x[i] - x_prev[i]));
            }
            if (max_norm < EPS) {
                std::cout << "Needed iterations: " << iter << std::endl;
                return x;
            }
        }
    }
};