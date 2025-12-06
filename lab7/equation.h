#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <utility>

#include "matrix.h"

enum SolverType { Exact, Seidel, Leibmann, SOR };

double norm_inf(const Matrix<double> &A);
double norm_inf_vec(const Matrix<double> &v);

struct EquationData {
    std::function<double(double)> phi0;
    std::function<double(double)> phi1;
    std::function<double(double)> phi2;
    std::function<double(double)> phi3;
    std::function<double(double, double)> exact;
    double l;
};

class Elliptic {
public:
    EquationData eq_data;
    std::function<std::string(int, Matrix<double> &, Matrix<double> &, double)> solve_func;
    double _h;
    std::ofstream fout;

public:
    Elliptic(const EquationData &data, const SolverType &solver);
    std::string solve(int N, double eps);

private:
    std::string leibmann_solve(int N, Matrix<double> &A, Matrix<double> &b, double eps);
    std::string seidel_solve(int N, Matrix<double> &A, Matrix<double> &b, double eps);
    std::string sor_solve(int N, Matrix<double> &A, Matrix<double> &b, double eps, double omega);
    std::string exact(int N);

    std::pair<Matrix<double>, Matrix<double>> get_equation_system(int N);
    std::pair<Matrix<double>, Matrix<double>> get_equivalent_system(Matrix<double> A, Matrix<double> b);
    Matrix<double> vector_to_matrix(const Matrix<double> &vec, int N);
};

inline Elliptic::Elliptic(const EquationData &data, const SolverType &type) : eq_data(data) {
    switch (type) {
        case Exact:
            solve_func = [this](int N, Matrix<double> &A, Matrix<double> &b, double eps) { return this->exact(N); };
            fout.open("exact.txt");
            break;
        case Seidel:
            solve_func = [this](int N, Matrix<double> &A, Matrix<double> &b, double eps) {
                return this->seidel_solve(N, A, b, eps);
            };
            fout.open("result.txt");
            std::cout << "seidel" << std::endl;
            break;
        case Leibmann:
            solve_func = [this](int N, Matrix<double> &A, Matrix<double> &b, double eps) {
                return this->leibmann_solve(N, A, b, eps);
            };
            fout.open("result.txt");
            std::cout << "Leibmann" << std::endl;
            break;
        case SOR:
            solve_func = [this](int N, Matrix<double> &A, Matrix<double> &b, double eps) {
                return this->sor_solve(N, A, b, eps, 1.8);
            };
            fout.open("result.txt");
            std::cout << "SOR" << std::endl;
            break;
        default:
            break;
    }
}

inline std::string Elliptic::exact(int N) {
    int sz = N - 1;

    for (int i = 0; i < sz; ++i) {
        for (int j = 0; j < sz; ++j) {
            double x = (i + 1) * _h;
            double y = (j + 1) * _h;
            fout << eq_data.exact(x, y) << " ";
        }
        fout << "\n";
    }
    fout << std::flush;
    return "exact";
}

inline std::string Elliptic::solve(int N, double eps) {
    _h = eq_data.l / N;
    auto [A, b] = get_equation_system(N);
    fout << N - 1 << ' ' << N - 1 << '\n';
    return solve_func(N, A, b, eps);
}

inline std::string Elliptic::seidel_solve(int N, Matrix<double> &A, Matrix<double> &b, double eps) {
    int n = A.getRows();
    int sz = N - 1;

    Matrix<double> alpha(n, n), beta(n, 1);
    std::tie(alpha, beta) = get_equivalent_system(A, b);

    double alpha_norm = 0.0;
    for (int i = 0; i < n; ++i) {
        double row_sum = 0.0;
        for (int j = 0; j < n; ++j) {
            row_sum += std::abs(alpha.at(i, j));
        }
        if (row_sum > alpha_norm) {
            alpha_norm = row_sum;
        }
    }

    double convergence_factor = 1.0;
    if (alpha_norm < 1.0) {
        convergence_factor = 1.0 - alpha_norm;
    } else {
        std::cerr << "ALPHA NORM IS MORE THEN ONE: " << alpha_norm << std::endl;
    }

    Matrix<double> x = beta;  // Начальное приближение
    Matrix<double> next_x(n, 1, 0.0);

    int iter = 0;
    int max_iter = 1000000;  // Увеличено для N >= 50

    while (true) {
        iter++;
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < i; ++j) sum += alpha.at(i, j) * next_x.at(j, 0);
            for (int j = i; j < n; ++j) sum += alpha.at(i, j) * x.at(j, 0);
            next_x.at(i, 0) = beta.at(i, 0) + sum;
        }

        double max_diff = 0.0;
        for (int i = 0; i < n; ++i) {
            double diff = std::abs(next_x.at(i, 0) - x.at(i, 0));
            if (diff > max_diff) max_diff = diff;
        }

        if (max_diff < eps * convergence_factor) {
            std::cout << "Iterative error: " << max_diff << std::endl;
            break;
        }

        if (iter > max_iter) {
            std::cout << "Max iterations reached! Alpha norm: " << alpha_norm << std::endl;
            break;
        }

        x = next_x;
    }

    Matrix<double> u = vector_to_matrix(x, sz);
    for (int i = 0; i < sz; ++i) {
        for (int j = 0; j < sz; ++j) {
            fout << u.at(i, j) << " ";
        }
        fout << "\n";
    }
    fout << std::flush;

    return "seidel iterations: " + std::to_string(iter);
}

inline std::pair<Matrix<double>, Matrix<double>> Elliptic::get_equation_system(int N) {
    int sz = N - 1;
    int size = sz * sz;
    double h_sq = _h * _h;

    Matrix<double> A(size, size, 0.0);
    Matrix<double> b(size, 1, 0.0);

    for (int i = 0; i < sz; ++i) {
        for (int j = 0; j < sz; ++j) {
            double x = (i + 1) * _h;
            double y = (j + 1) * _h;

            int row = i * sz + j;

            A.at(row, row) = -4.0 + 2.0 * h_sq;

            if (i + 1 == sz) {  // x = L
                b.at(row, 0) -= eq_data.phi1(y);
            } else {
                A.at(row, (i + 1) * sz + j) = 1.0;
            }

            if (i - 1 == -1) {  // x = 0
                b.at(row, 0) -= eq_data.phi0(y);
            } else {
                A.at(row, (i - 1) * sz + j) = 1.0;
            }

            if (j + 1 == sz) {  // y = L
                b.at(row, 0) -= eq_data.phi3(x);
            } else {
                A.at(row, i * sz + (j + 1)) = 1.0;
            }

            if (j - 1 == -1) {  // y = 0
                b.at(row, 0) -= eq_data.phi2(x);
            } else {
                A.at(row, i * sz + (j - 1)) = 1.0;
            }
        }
    }
    return {A, b};
}

inline std::string Elliptic::leibmann_solve(int N, Matrix<double> &A, Matrix<double> &b, double eps) {
    int n = A.getRows();
    int sz = N - 1;

    Matrix<double> alpha(n, n), beta(n, 1);
    std::tie(alpha, beta) = get_equivalent_system(A, b);

    double alpha_norm = 0.0;
    for (int i = 0; i < n; ++i) {
        double row_sum = 0.0;
        for (int j = 0; j < n; ++j) {
            row_sum += std::abs(alpha.at(i, j));
        }
        if (row_sum > alpha_norm) {
            alpha_norm = row_sum;
        }
    }

    double convergence_factor = 1.0;
    if (alpha_norm < 1.0) {
        convergence_factor = 1.0 - alpha_norm;
    } else {
        std::cerr << "ALPHA NORM IS MORE THEN ONE: " << alpha_norm << std::endl;
    }

    Matrix<double> x = beta;
    Matrix<double> next_x(n, 1, 0.0);
    int iter = 0;
    while (true) {
        iter++;
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                sum += alpha.at(i, j) * x.at(j, 0);
            }
            next_x.at(i, 0) = beta.at(i, 0) + sum;
        }
        double max_diff = 0.0;
        for (int i = 0; i < n; ++i) {
            double d = std::abs(next_x.at(i, 0) - x.at(i, 0));
            if (d > max_diff) max_diff = d;
        }

        if (max_diff < eps * convergence_factor) {
            std::cout << "Iteration error: " << max_diff << std::endl;
            break;
        }

        if (iter > 1000000) {
            std::cout << "Leibmann max iter reached\n";
            break;
        }

        x = next_x;
    }

    Matrix<double> u = vector_to_matrix(x, sz);
    for (int i = 0; i < sz; ++i) {
        for (int j = 0; j < sz; ++j) {
            fout << u.at(i, j) << " ";
        }
        fout << "\n";
    }
    fout << std::flush;

    return "leibmann iterations: " + std::to_string(iter);
}

inline std::pair<Matrix<double>, Matrix<double>> Elliptic::get_equivalent_system(Matrix<double> A, Matrix<double> b) {
    int n = A.getRows();

    Matrix<double> alpha(n, n, 0.0);
    Matrix<double> beta(n, 1, 0.0);

    for (int i = 0; i < n; ++i) {
        if (A.at(i, i) == 0.0) {
            bool swapped = false;
            for (int j = i + 1; j < n; ++j) {
                if (A.at(j, i) != 0.0) {
                    for (int k = 0; k < n; ++k) {
                        std::swap(A.at(i, k), A.at(j, k));
                    }
                    std::swap(b.at(i, 0), b.at(j, 0));
                    swapped = true;
                    break;
                }
            }
            if (!swapped) {
                continue;  // аналог Python: просто пропустить эту строку
            }
        }
        beta.at(i, 0) = b.at(i, 0) / A.at(i, i);
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                alpha.at(i, j) = 0.0;
            } else {
                alpha.at(i, j) = -A.at(i, j) / A.at(i, i);
            }
        }
    }

    return {alpha, beta};
}

inline std::string Elliptic::sor_solve(int N, Matrix<double> &A, Matrix<double> &b, double eps, double omega) {
    int n = A.getRows();
    int sz = N - 1;

    Matrix<double> alpha(n, n), beta(n, 1);
    std::tie(alpha, beta) = get_equivalent_system(A, b);

    Matrix<double> x = beta;
    Matrix<double> next_x(n, 1, 0.0);
    next_x = x;

    int iter = 0;
    int max_iter = 1000000;

    while (true) {
        iter++;
        double max_diff = 0.0;

        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < i; ++j) sum += alpha.at(i, j) * next_x.at(j, 0);
            for (int j = i; j < n; ++j) sum += alpha.at(i, j) * x.at(j, 0);
            double u_GS = beta.at(i, 0) + sum;
            double x_new = omega * u_GS + (1.0 - omega) * x.at(i, 0);
            double diff = std::abs(x_new - x.at(i, 0));
            if (diff > max_diff) max_diff = diff;
            next_x.at(i, 0) = x_new;
        }
        if (max_diff < eps) {
            std::cout << "Iteration error: " << max_diff << std::endl;
            break;
        }

        if (iter > max_iter) {
            std::cout << "SOR max iterations reached! Omega: " << omega << std::endl;
            break;
        }
        x = next_x;
    }

    Matrix<double> u = vector_to_matrix(x, sz);
    for (int i = 0; i < sz; ++i) {
        for (int j = 0; j < sz; ++j) {
            fout << u.at(i, j) << " ";
        }
        fout << "\n";
    }
    fout << std::flush;

    return "sor, iterations: " + std::to_string(iter);
}

inline Matrix<double> Elliptic::vector_to_matrix(const Matrix<double> &vec, int sz) {
    Matrix<double> u(sz, sz, 0.0);

    for (int i = 0; i < sz; ++i) {
        for (int j = 0; j < sz; ++j) {
            u.at(i, j) = vec.at(i * sz + j, 0);
        }
    }

    return u;
}

inline double norm_inf(const Matrix<double> &A) {
    int rows = A.getRows();
    int cols = A.getCols();

    double max_sum = 0.0;

    for (int i = 0; i < rows; ++i) {
        double row_sum = 0.0;
        for (int j = 0; j < cols; ++j) {
            row_sum += std::abs(A.at(i, j));
        }
        if (row_sum > max_sum) max_sum = row_sum;
    }

    return max_sum;
}

inline double norm_inf_vec(const Matrix<double> &v) {
    int n = v.getRows();

    double max_val = 0.0;

    for (int i = 0; i < n; ++i) {
        double val = std::abs(v.at(i, 0));
        if (val > max_val) max_val = val;
    }

    return max_val;
}