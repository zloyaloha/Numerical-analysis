#include <iostream>
#include "matrix.h"
#include <vector>
#include <stdexcept>
#include <tgmath.h>
#include <math.h>

const double EPSILON = 0.001;

struct MatrixElem {
    int row;
    int col;
    long double value;
};

class Rotating {
public:
    Rotating(long double epsilon) : eps(epsilon) {}

    std::pair<std::vector<long double>, Matrix<long double>> solve(Matrix<long double> matrix);
private:
    MatrixElem findMaxOverDiagonal(const Matrix<long double>& matrix);
    long double calculateNonDiagonalSummary(const Matrix<long double>& matrix);
    long double eps;
};

std::pair<std::vector<long double>, Matrix<long double>> Rotating::solve(Matrix<long double> matrix)
{
    std::vector<Matrix<long double>> Us;
    int k = 0;

    while(calculateNonDiagonalSummary(matrix) > eps) {
        MatrixElem max_elem = findMaxOverDiagonal(matrix);
        long double phi, cos_phi, sin_phi;
        if (std::fabs(matrix.at(max_elem.col, max_elem.col) - matrix.at(max_elem.row, max_elem.row)) < eps) {
            cos_phi = 1 / sqrtl(2);
            sin_phi = 1 / sqrtl(2);
        } else {
            long double theta = (matrix.at(max_elem.row, max_elem.row) - matrix.at(max_elem.col, max_elem.col)) / (2.0 * matrix.at(max_elem.col, max_elem.row));
            if (std::fabs(theta) < 1e-6) {
                theta = copysign(1e-6, theta);  // Минимальный угол, чтобы не зависнуть
            }
            long double t = copysign(1.0 / (fabs(theta) + sqrtl(1.0 + theta * theta)), theta);
            cos_phi = 1 / sqrtl(1 + t * t);
            sin_phi = t * cos_phi;
        }

        Matrix<long double> U(matrix.getCols(), matrix.getRows(), 0);
        for (int i = 0; i < matrix.getRows(); ++i) {
            U.at(i, i) = 1;
        }

        Us.push_back(U);

        Us[k].at(max_elem.row, max_elem.row) = cos_phi; Us[k].at(max_elem.row, max_elem.col) = -sin_phi;
        Us[k].at(max_elem.col, max_elem.row) = sin_phi; Us[k].at(max_elem.col, max_elem.col) = cos_phi;
        matrix = ~Us[k] * matrix * Us[k];

        if (k > 100 * matrix.getRows()) {
            break;
        }
        ++k;
    }

    if (calculateNonDiagonalSummary(matrix) < eps) {
        std::cout << "Result for accuracy = " << eps << " by " << k << " iterations" << std::endl;
    } else {
        std::cout << "Can't reach accuracy = " << eps << " by " << k << " iterations" << std::endl;
    }

    std::cout << "eigenvalues \n\n";
    std::vector<long double> eigenvalues;
    for (int i = 0; i < matrix.getCols(); ++i) {
        eigenvalues.push_back(matrix.at(i, i));
        std::cout << eigenvalues[i] << std::endl;
    }
    Matrix<long double> eigenvectors = Us[0] * Us[1];

    for (int i = 2; i < Us.size(); ++i) {
        eigenvectors = eigenvectors * Us[i];
    }

    std::cout << "\neigenvectors \n\n";
    eigenvectors.print();

    return std::make_pair(eigenvalues, eigenvectors);
}

MatrixElem Rotating::findMaxOverDiagonal(const Matrix<long double>& matrix) {
    MatrixElem max_elem;
    max_elem.value = 0;
    for (int i = 0; i < matrix.getCols(); ++i) {
        for (int j = i + 1; j < matrix.getRows(); ++j) {
            if (max_elem.value < std::fabs(matrix.at(i, j))) {
                max_elem.value = std::fabs(matrix.at(i, j));
                max_elem.row = i;
                max_elem.col = j;
            }
        }
    }
    return max_elem;
}

long double Rotating::calculateNonDiagonalSummary(const Matrix<long double>& matrix) {
    long double summary = 0;
    for (int i = 0; i < matrix.getRows(); ++i) {
        for (int j = i + 1; j < matrix.getRows(); ++j) {
            summary += std::pow(matrix.at(i, j), 2);
        }
    }
    return std::pow(summary, 0.5);
}

