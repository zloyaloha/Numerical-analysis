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

double characteristicRoots(double a, double b, double c, double d, double alpha) {
    double trace = a + d;
    double det = a * d - b * c;

    double discriminant = trace * trace - 4 * det;

    if (discriminant < 0) {
        double realPart = trace / 2;
        double imagPart = sqrt(-discriminant) / 2;
        std::cout << "Корни характеристического уравнения (комплексные):\n";
        std::cout << "λ1 = " << realPart << " + " << imagPart << "i\n";
        std::cout << "λ2 = " << realPart << " - " << imagPart << "i\n";
        std::cout << "Корень характеристического уравнения (не комплексный):\n";
        std::cout << "λ2 = " << alpha << std::endl;
    }
    return discriminant;
}

class QR {
private:
    std::pair<Matrix<double>, Matrix<double>> householder(Matrix<double> A);
    bool isConverged(const Matrix<double>& A, double eps);
public:
    QR() = default;
    void decompose(Matrix<double> A, const double& eps);
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
        Matrix<double> x(n - k, 1, 0.0);
        for (int i = k; i < n; ++i)
            x.at(i - k, 0) = R.at(i, k);

        Matrix<double> e(n - k, 1, 0.0);
        e.at(0, 0) = 1.0;

        double alpha = norm(x);
        if (x.at(0, 0) >= 0)
            alpha = -alpha;

        Matrix<double> v = x - e * alpha;

        Matrix<double> H_sub = Matrix<double>::identity(n - k) - (v * ~v) * (2.0 / ((~v * v).at(0, 0)));

        Matrix<double> H = Matrix<double>::identity(n);
        for (int i = k; i < n; ++i)
            for (int j = k; j < n; ++j)
                H.at(i, j) = H_sub.at(i - k, j - k);

        R = H * R;
        Q = Q * H;
    }

    return {Q, R};
}

void QR::decompose(Matrix<double> A, const double& eps) {
    int n = A.getCols();
    int iter = 0;
    while (true) {
        ++iter;
        auto [Q, R] = householder(A);
        A = (R * Q);
        if (isConverged(A, eps)) {
            std::cout << "EIGENVALUES!\n";
            for (int i = 0; i < A.getCols(); ++i) {
                std::cout << A.at(i, i) << '\n';
            }
            return;
        }
        if (std::sqrt(A.at(2, 0) * A.at(2, 0)) < eps && (std::sqrt(A.at(1, 0) * A.at(1, 0)) > eps || std::sqrt(A.at(2, 1) * A.at(2, 1)) > eps)) {
            double discr;
            if (std::sqrt(A.at(1, 0) * A.at(1, 0)) > eps) {
                discr = characteristicRoots(A.at(0, 0), A.at(0, 1), A.at(1, 0), A.at(0, 1), A.at(2, 2));
            }
            if (std::sqrt(A.at(2, 1) * A.at(2, 1)) > eps) {
                discr = characteristicRoots(A.at(1, 1), A.at(1, 2), A.at(2, 1), A.at(2, 2), A.at(0, 0));
            }
            if (discr < 0) {
                return;
            }
        }
        if (iter == MAX_ITER) {
            std::cerr << "Не удалось достичь желаемой точности за " << iter << " итераций\n";
            return;
        }
    }
}