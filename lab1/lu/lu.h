#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include "matrix.h"
#include <iomanip>
#include <stdexcept>

class LU {

public:
    LU() = default;

    template <class T>
    typename std::enable_if<std::is_arithmetic<T>::value, std::pair<Matrix<double>, Matrix<double>>>::type
    decompose(Matrix<T>& A, bool info = 0)
    {
        Matrix<double> L(A.getRows(), A.getRows(), 0.), U(A.getRows(), A.getRows(), 0.);
        std::vector<int> P(A.getRows());
        for (int i = 0; i < A.getRows(); ++i) {
            P[i] = i;
        }
        int permutatuions = 0;

        for (int k = 0; k < A.getRows(); ++k) {
            double max_elem = abs(A.at(k, k));
            int max_ind = k;
            for (int i = k + 1; i < A.getRows(); ++i) {
                if (abs(A.at(i, k)) > max_elem) {
                    max_elem = abs(A.at(i, k));
                    max_ind = i;
                }
            }

            if (max_elem == 0) {
                throw std::invalid_argument("Broken matrix!");
            }

            if (max_ind != k) {
                ++permutatuions;
                A.swapRows(k, max_ind);
                std::swap(P[k], P[max_ind]);
            }
        }

        for (int k = 0; k < A.getRows(); ++k) {
            U.at(k, k) = A.at(k, k);
            L.at(k, k) = 1;
            for (int i = k + 1; i < A.getRows(); ++i) {
                L.at(i, k) = A.at(i, k) / U.at(k, k);
                U.at(k, i) = A.at(k, i);
            }

            for (int i = k + 1; i < A.getRows(); ++i) {
                for (int j = k + 1; j < A.getRows(); ++j) {
                    A.at(i, j) -= L.at(i, k) * U.at(k, j);
                }
            }
        }
        if (info) {
            std::cout << "determinant: " << determinant(L, U) * (permutatuions % 2 == 0 ? 1 : -1) << '\n';
            std::cout << "inverted\n";
            Matrix<double> AInv = invert(L, U, P);
            AInv.print();
        }
        return std::make_pair(L, U);
    }

    template <class T>
    typename std::enable_if<std::is_arithmetic<T>::value, std::vector<double>>::type
    solve(Matrix<T> A)
    {
        auto [L, U] = decompose(A);

        std::vector<double> x;
        for (int i = 0; i < A.getRows(); ++i) {
            x.push_back(A.at(i, A.getCols() - 1));
        }

        return backward_substitution(U, forward_substitution(L, x));
    }

    private:
    std::vector<double> forward_substitution(const Matrix<double>& L, const std::vector<double>& b)
    {
        size_t n = L.getRows();
        std::vector<double> y(n, 0);

        for (size_t i = 0; i < n; ++i) {
            y[i] = b[i];
            for (size_t j = 0; j < i; ++j) {
                y[i] -= L.at(i, j) * y[j];
            }
        }
        return y;
    }

    std::vector<double> backward_substitution(const Matrix<double>& U, const std::vector<double>& y)
    {
        size_t n = U.getRows();
        std::vector<double> x(n, 0);

        for (int i = n - 1; i >= 0; --i) {
            x[i] = y[i];
            for (size_t j = i + 1; j < n; ++j) {
                x[i] -= U.at(i, j) * x[j];
            }
            x[i] /= U.at(i, i);
        }

        return x;
    }

    double determinant(const Matrix<double>& L, const Matrix<double>& U)
    {
        double det = 1;
        for (int i = 0; i < U.getRows(); ++i) {
            det *= (U.at(i, i) * L.at(i, i));
        }
        return det;
    }

    Matrix<double> invert(const Matrix<double>& L, const Matrix<double>& U, const std::vector<int> P)
    {
        Matrix<double> invA(L.getRows(), L.getRows());

        for (int i = 0; i < L.getRows(); i++) {
            std::vector<double> e(L.getRows(), 0.0);
            e[i] = 1.0;

            std::vector<double> Y = forward_substitution(L, e);
            std::vector<double> X = backward_substitution(U, Y);

            for (int j = 0; j < L.getRows(); j++)
                invA.at(j, P[i]) = X[j];
        }
        return invA;
    }
};