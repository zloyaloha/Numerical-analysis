#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "tridiagonal-matrix.h"

enum SolverType { Exact, AlterDirection, FractSteps };

struct EquationData {
    std::function<double(double, double, double)> f;
    std::function<double(double, double)> phi0;
    std::function<double(double, double)> phi1;
    std::function<double(double, double)> phi2;
    std::function<double(double, double)> phi3;
    std::function<double(double, double)> psi;
    std::function<double(double, double, double)> exact;
    double l1, l2;
};

class Parabolic2D {
public:
    EquationData eq_data;
    std::function<std::string(int, int, int, double)> solve_func;
    double _h1, _h2, _tau;

    int _n_snapshots, _save_step;
    std::ofstream fout;

public:
    Parabolic2D(const EquationData& data, const SolverType& solver, int n_snapshots = 5);
    ~Parabolic2D() {
        if (fout.is_open()) fout.close();
    }
    std::string solve(int N1, int N2, double sigma, double T);

private:
    std::string fract_steps(int N1, int N2, int K, double T);
    std::string alter_direction(int N1, int N2, int K, double T);
    std::string exact(int N1, int N2, int K, double T);

    void save_snapshot(const std::vector<std::vector<double>>& u, double t);
    double calc_max_error(const std::vector<std::vector<double>>& u, double t);
};

inline Parabolic2D::Parabolic2D(const EquationData& data, const SolverType& type, int n_snapshots)
    : eq_data(data), _n_snapshots(n_snapshots) {
    switch (type) {
        case FractSteps:
            solve_func = [this](int N1, int N2, int K, double T) { return this->fract_steps(N1, N2, K, T); };
            fout.open("result_fract_steps.txt");
            std::cout << "Solver: Fract Steps" << std::endl;
            break;
        case AlterDirection:
            solve_func = [this](int N1, int N2, int K, double T) { return this->alter_direction(N1, N2, K, T); };
            fout.open("result_adi.txt");
            std::cout << "Solver: ADI" << std::endl;
            break;
        case Exact:
            solve_func = [this](int N1, int N2, int K, double T) { return this->exact(N1, N2, K, T); };
            fout.open("result_exact.txt");
            std::cout << "Solver: Exact" << std::endl;
            break;
        default:
            solve_func = [](int, int, int, double) { return "Error: Solver not found."; };
            break;
    }
}

inline std::string Parabolic2D::solve(int N1, int N2, double sigma, double T) {
    _h1 = eq_data.l1 / N1;
    _h2 = eq_data.l2 / N2;

    double h_min = std::min(_h1, _h2);

    _tau = sigma * h_min * h_min;

    int K = std::ceil(T / _tau);
    _tau = T / K;

    _save_step = std::max(1, K / _n_snapshots);

    std::cout << "Grid: " << N1 << "x" << N2 << ", Steps: " << K << ", Save every: " << _save_step << " steps." << std::endl;

    return solve_func(N1, N2, K, T);
}

inline std::string Parabolic2D::exact(int N1, int N2, int K, double T) {
    for (int i = 0; i <= N1; ++i) {
        double x = i * _h1;
        for (int j = 0; j <= N2; ++j) {
            double y = j * _h2;
            double val = eq_data.exact(x, y, T);
            fout << std::fixed << std::setprecision(10) << val << (j == N2 ? "" : " ");
        }
        fout << "\n";
    }
    return "Exact solution written.";
}

inline std::string Parabolic2D::fract_steps(int N1, int N2, int K, double T) {
    std::vector<std::vector<double>> u_curr(N1 + 1, std::vector<double>(N2 + 1));
    std::vector<std::vector<double>> u_half(N1 + 1, std::vector<double>(N2 + 1));
    std::vector<std::vector<double>> u_next(N1 + 1, std::vector<double>(N2 + 1));

    std::ofstream error_file("fract_steps_error.txt");

    for (int i = 0; i <= N1; ++i) {
        for (int j = 0; j <= N2; ++j) {
            u_curr[i][j] = eq_data.psi(i * _h1, j * _h2);
        }
    }
    error_file << calc_max_error(u_curr, 0) << ' ';

    double r_x = _tau / (_h1 * _h1);
    double r_y = _tau / (_h2 * _h2);

    for (int k = 0; k < K; ++k) {
        double t_curr = k * _tau;
        double t_next = (k + 1) * _tau;
        double t_half = t_curr + 0.5 * _tau;

        for (int j = 1; j < N2; ++j) {
            std::vector<double> a(N1 + 1, 0.0), c(N1 + 1, 0.0), b(N1 + 1, 0.0), f_rhs(N1 + 1, 0.0);
            std::vector<double> x_res(N1 + 1);

            double y = j * _h2;

            c[0] = 1.0;
            b[0] = 0.0;
            f_rhs[0] = eq_data.phi0(y, t_half);

            for (int i = 1; i < N1; ++i) {
                double x = i * _h1;
                a[i] = -r_x;
                b[i] = -r_x;
                c[i] = 1.0 + 2.0 * r_x;

                f_rhs[i] = u_curr[i][j] + (_tau / 2.0) * eq_data.f(x, y, t_curr);
            }
            a[N1] = 0.0;
            c[N1] = 1.0;
            f_rhs[N1] = eq_data.phi1(y, t_half);

            x_res = tma(a, c, b, f_rhs);

            for (int i = 0; i <= N1; ++i) u_half[i][j] = x_res[i];
        }

        for (int i = 0; i <= N1; ++i) {
            u_half[i][0] = eq_data.phi2(i * _h1, t_half);
            u_half[i][N2] = eq_data.phi3(i * _h1, t_half);
        }

        for (int i = 1; i < N1; ++i) {
            std::vector<double> a(N2 + 1, 0.0), c(N2 + 1, 0.0), b(N2 + 1, 0.0), f_rhs(N2 + 1, 0.0);
            std::vector<double> y_res(N2 + 1);

            double x = i * _h1;

            c[0] = 1.0;
            b[0] = 0.0;
            f_rhs[0] = eq_data.phi2(x, t_next);

            for (int j = 1; j < N2; ++j) {
                double y = j * _h2;
                a[j] = -r_y;
                b[j] = -r_y;
                c[j] = 1.0 + 2.0 * r_y;

                f_rhs[j] = u_half[i][j] + (_tau / 2.0) * eq_data.f(x, y, t_curr);
            }
            a[N2] = 0.0;
            c[N2] = 1.0;
            f_rhs[N2] = eq_data.phi3(x, t_next);

            y_res = tma(a, c, b, f_rhs);

            for (int j = 0; j <= N2; ++j) u_next[i][j] = y_res[j];
        }

        for (int j = 0; j <= N2; ++j) {
            u_next[0][j] = eq_data.phi0(j * _h2, t_next);
            u_next[N1][j] = eq_data.phi1(j * _h2, t_next);
        }
        for (int i = 0; i <= N1; ++i) {
            u_next[i][0] = eq_data.phi2(i * _h1, t_next);
            u_next[i][N2] = eq_data.phi3(i * _h1, t_next);
        }

        u_curr = u_next;

        error_file << calc_max_error(u_curr, t_next) << ' ';
        if ((k + 1) % _save_step == 0 || (k + 1) == K) {
            save_snapshot(u_curr, t_next);
        }
    }

    return "Splitting method solution";
}

inline std::string Parabolic2D::alter_direction(int N1, int N2, int K, double T) {
    std::vector<std::vector<double>> u_curr(N1 + 1, std::vector<double>(N2 + 1));
    std::vector<std::vector<double>> u_half(N1 + 1, std::vector<double>(N2 + 1));  // слой k + 1/2
    std::vector<std::vector<double>> u_next(N1 + 1, std::vector<double>(N2 + 1));
    std::ofstream error_file("adi_error.txt");

    for (int i = 0; i <= N1; ++i)
        for (int j = 0; j <= N2; ++j) u_curr[i][j] = eq_data.psi(i * _h1, j * _h2);
    error_file << calc_max_error(u_curr, 0) << ' ';

    double r_x = _tau / (2.0 * _h1 * _h1);
    double r_y = _tau / (2.0 * _h2 * _h2);

    int n_inner_x = N1 - 1;
    std::vector<double> ax(n_inner_x, -r_x);
    std::vector<double> bx(n_inner_x, 1.0 + 2.0 * r_x);
    std::vector<double> cx(n_inner_x, -r_x);
    std::vector<double> dx(n_inner_x);

    int n_inner_y = N2 - 1;
    std::vector<double> ay(n_inner_y, -r_y);
    std::vector<double> by(n_inner_y, 1.0 + 2.0 * r_y);
    std::vector<double> cy(n_inner_y, -r_y);
    std::vector<double> dy(n_inner_y);

    for (int k = 0; k < K; ++k) {
        double t_curr = k * _tau;
        double t_half = (k + 0.5) * _tau;
        double t_next = (k + 1) * _tau;

        for (int j = 0; j <= N2; ++j) {
            u_half[0][j] = eq_data.phi0(j * _h2, t_half);
            u_half[N1][j] = eq_data.phi1(j * _h2, t_half);
        }
        for (int i = 0; i <= N1; ++i) {
            u_half[i][0] = eq_data.phi2(i * _h1, t_half);
            u_half[i][N2] = eq_data.phi3(i * _h1, t_half);
        }

        for (int j = 1; j < N2; ++j) {
            for (int i = 1; i < N1; ++i) {
                double term_y = r_y * (u_curr[i][j + 1] - 2.0 * u_curr[i][j] + u_curr[i][j - 1]);
                double source = (_tau / 2.0) * eq_data.f(i * _h1, j * _h2, t_half);

                dx[i - 1] = u_curr[i][j] + term_y + source;
            }
            dx[0] -= (-r_x) * u_half[0][j];
            dx[n_inner_x - 1] -= (-r_x) * u_half[N1][j];

            std::vector<double> sol_x = tma(ax, bx, cx, dx);
            for (int i = 1; i < N1; ++i) u_half[i][j] = sol_x[i - 1];
        }

        for (int j = 0; j <= N2; ++j) {
            u_next[0][j] = eq_data.phi0(j * _h2, t_next);
            u_next[N1][j] = eq_data.phi1(j * _h2, t_next);
        }
        for (int i = 0; i <= N1; ++i) {
            u_next[i][0] = eq_data.phi2(i * _h1, t_next);
            u_next[i][N2] = eq_data.phi3(i * _h1, t_next);
        }

        for (int i = 1; i < N1; ++i) {
            for (int j = 1; j < N2; ++j) {
                double term_x = r_x * (u_half[i + 1][j] - 2.0 * u_half[i][j] + u_half[i - 1][j]);
                double source = (_tau / 2.0) * eq_data.f(i * _h1, j * _h2, t_half);

                dy[j - 1] = u_half[i][j] + term_x + source;
            }

            dy[0] -= (-r_y) * u_next[i][0];
            dy[n_inner_y - 1] -= (-r_y) * u_next[i][N2];

            std::vector<double> sol_y = tma(ay, by, cy, dy);
            for (int j = 1; j < N2; ++j) u_next[i][j] = sol_y[j - 1];
        }

        u_curr = u_next;
        error_file << calc_max_error(u_curr, t_next) << ' ';

        if ((k + 1) % _save_step == 0 || (k + 1) == K) {
            save_snapshot(u_curr, t_next);
        }
    }
    return "ADI solution";
}

inline void Parabolic2D::save_snapshot(const std::vector<std::vector<double>>& u, double t) {
    int N1 = u.size() - 1;
    int N2 = u[0].size() - 1;

    fout << "TIME " << t << "\n";

    for (int i = 0; i <= N1; ++i) {
        for (int j = 0; j <= N2; ++j) {
            fout << (i * _h1) << " " << (j * _h2) << " " << u[i][j] << "\n";
        }
    }
    fout << "\n";
}

inline double Parabolic2D::calc_max_error(const std::vector<std::vector<double>>& u, double t) {
    std::vector<double> errors;
    double max_err = 0.0;
    int N1 = u.size() - 1;
    int N2 = u[0].size() - 1;

    for (int i = 0; i <= N1; ++i) {
        for (int j = 0; j <= N2; ++j) {
            double exact_val = eq_data.exact(i * _h1, j * _h2, t);
            double diff = std::abs(u[i][j] - exact_val);
            max_err = std::fmax(diff, max_err);
        }
    }
    return max_err;
}