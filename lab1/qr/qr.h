#include "matrix.h"
#include <cmath>

const double EPS = 1e-10;
const int MAX_ITER = 100;

double norm(const Matrix<double>& v) {
    double sum = 0.0;
    for (int i = 0; i < v.getRows(); ++i) {
        sum += v.at(i, 0) * v.at(i, 0);
    }
    return std::sqrt(sum);
}


class QR {
private:
    std::pair<Matrix<double>, Matrix<double>> householder(Matrix<double> A);
    bool isConverged(const Matrix<double>& A, double eps);
public:
    QR() = default;
    Matrix<double>  decompose(Matrix<double> A, const double& eps);
};

bool QR::isConverged(const Matrix<double>& A, double eps) {
    int n = A.getCols();
    double sum = 0.0;
    for (int i = 1; i < n; ++i)
        for (int j = 0; j < i; ++j)
            sum += std::pow(A.at(i, j), 2);
    return sqrt(sum) < eps;
}

std::pair<Matrix<double>, Matrix<double>> QR::householder(Matrix<double> A)
{
    int n = A.getRows();
    Matrix<double> R = A;
    Matrix<double> Q = Matrix<double>::identity(n);

    for (int k = 0; k < n - 1; ++k) {
        // Подвектор x (только нижняя часть столбца)
        Matrix<double> x(n - k, 1, 0.0);
        for (int i = k; i < n; ++i)
            x.at(i - k, 0) = R.at(i, k);

        // Вектор e
        Matrix<double> e(n - k, 1, 0.0);
        e.at(0, 0) = 1.0;

        // Вычисление alpha
        double alpha = norm(x);
        if (x.at(0, 0) >= 0)
            alpha = -alpha;

        // Вектор v = x - alpha * e
        Matrix<double> v = x - e * alpha;

        // Нормализация v
        double v_norm = norm(v);
        if (v_norm > 1e-12) {
            v = v * (1.0 / v_norm);
        }

        // Матрица Хаусхолдера H_sub = I - 2vv^T
        Matrix<double> H_sub = Matrix<double>::identity(n - k) - (v * ~v) * (2.0 / ((~v * v).at(0, 0)));

        // Встраиваем H_sub в полную H
        Matrix<double> H = Matrix<double>::identity(n);
        for (int i = k; i < n; ++i)
            for (int j = k; j < n; ++j)
                H.at(i, j) = H_sub.at(i - k, j - k);

        R = H * R;
        Q = Q * H;
    }

    return {Q, R};
}

Matrix<double> QR::decompose(Matrix<double> A, const double& eps) {
    int n = A.getCols();
    int iter = 0;
    while (true) {
        ++iter;
        auto [Q, R] = householder(A);
        A = (R * Q);
        if (isConverged(A, eps)) {
            return A;
        }
        if (iter == MAX_ITER) {
            std::cerr << "Не удалось достичь желаемой точности за " << iter << " итераций\n";
            return A;
        }

        A.print(); std::cout << std::endl;
    }
}