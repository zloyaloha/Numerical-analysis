#include <functional>
#include <fstream>
#include "tridiagonal-matrix.h"
#include <iostream>
#include <string>
#include <vector>

enum SolverType {
    Explicit,
    Implicit,
    Exact
};

struct EquationData {
    std::function<double(double, double)> f;
    std::function<double(double)> psi1;
    std::function<double(double)> psi2;
    std::function<double(double)> psi1_dir1;
    std::function<double(double)> psi1_dir2;
    std::function<double(double)> phi0;
    std::function<double(double)> phil;
    std::function<double(double, double)> exact_solution;
    double a, b, c, d, l;
    double alpha, beta, gamma, delta;
    std::string bound_type;
    std::string approximation;
};

class Hyperbolic {
public:
    EquationData eq_data;
    std::function<std::string(int,int)> solve_func;
    std::function<double(std::vector<std::vector<double>>&, double K, double T)> left_bound;
    std::function<double(std::vector<std::vector<double>>&, double K, double T)> right_bound;
    double _h, _tau, _sigma, _a, _omega;
    std::ofstream fout;
public:
    Hyperbolic(const EquationData& data, const SolverType& solver);
    std::string solve(int N, double sigma, double T, double a);
private:
    std::string explicit_solve(int N, int K);
    std::string implicit_solve(int N, int K);
    std::string exact(int N, int K);
};

inline Hyperbolic::Hyperbolic(const EquationData& data, const SolverType& type) : eq_data(data)
{
    switch (type)
    {
    case Explicit:
        solve_func = [this](int N, int K) { return this->explicit_solve(N, K); };
        fout.open("result.txt");
        std::cout << "explicit" << std::endl;
        break;
    case Implicit:
        solve_func = [this](int N, int K) { return this->implicit_solve(N, K); };
        fout.open("result.txt");
        std::cout << "impicit" << std::endl;
        break;
    case Exact:
        solve_func = [this](int N, int K) { return this->exact(N, K); };
        fout.open("exact.txt");
        break;
    default:
        break;
    }
}

inline std::string Hyperbolic::exact(int N, int K)
{
     for (int k = 0; k < K; ++k) {
        for (int j = 0; j < N; ++j) {
            double x = j * _h;
            double t = k * _tau;
            fout << eq_data.exact_solution(x, t) << " ";
        }
        fout << "\n";
    }

    fout << std::flush;
    return "exact";
}

inline std::string Hyperbolic::solve(int N, double sigma, double T, double a)
{
    _sigma = sigma;

    _h = eq_data.l / N;

    _tau = std::sqrt(_sigma) * _h / std::sqrt(a);
    _omega = _tau * _tau * eq_data.b / (2.0 * _h);


    int K = std::ceil(T / _tau);
    fout << K << ' ' << N << '\n';
    return solve_func(N, K);
}

inline std::string Hyperbolic::explicit_solve(int N, int K)
{
    std::vector<std::vector<double>> u(K, std::vector<double>(N));

    for (int j = 0; j < N; ++j) {
        double x = j * _h;
        u[0][j] = eq_data.psi1(x);

        double O = _tau * _tau / 2;
        if (eq_data.approximation == "p1") {
            u[1][j] = eq_data.psi1(x) + _tau * eq_data.psi2(x) + eq_data.psi1_dir2(x) * O;
        } else if (eq_data.approximation == "p2") {
            u[1][j] = eq_data.psi1(x) + _tau * eq_data.psi2(x) +
                eq_data.a * O * eq_data.psi1_dir2(x) +
                eq_data.b * O * eq_data.psi1_dir1(x) +
                eq_data.c * O * eq_data.psi1(x) +
                eq_data.d * O * eq_data.psi2(x) +
                O * eq_data.f(x, 0);
        }
    }

    double tau = _tau;
    u[1][0] = eq_data.phi0(tau);
    u[1][N - 1] = eq_data.phil(tau);

    for (int k = 2; k < K; ++k) {
        double t = (k - 1) * _tau;

        double tau_sq = _tau * _tau;
        double d_tau_half = eq_data.d * _tau / 2.0;

        double sigma_a_coeff = _sigma;
        double sigma_b_coeff = eq_data.b * tau_sq / (2.0 * _h);

        for (int j = 1; j < N - 1; ++j) {
            double x = j * _h;

            double numerator =
                (2.0 + eq_data.c * tau_sq - 2.0 * sigma_a_coeff) * u[k - 1][j] +
                -(1.0 - d_tau_half) * u[k - 2][j] +
                sigma_a_coeff * (u[k - 1][j + 1] + u[k - 1][j - 1]) +
                sigma_b_coeff * (u[k - 1][j + 1] - u[k - 1][j - 1]) +
                tau_sq * eq_data.f(x, t);
                u[k][j] = numerator / (1.0 + d_tau_half);
        }
        u[k][0] = eq_data.phi0((k) * _tau);
        u[k][N - 1] = eq_data.phil((k) * _tau);
    }

    for (int k = 0; k < K; ++k) {
        for (int j = 0; j < N; ++j) {
            fout << u[k][j] << " ";
        }
        fout << "\n";
    }
    fout << std::flush;
    fout.close();

    return "explicit";
}

inline std::string Hyperbolic::implicit_solve(int N, int K)
{
    std::vector<std::vector<double>> u(K, std::vector<double>(N, 0));

    for (int j = 0; j < N; ++j) {
        double x = j * _h;
        u[0][j] = eq_data.psi1(x);

        double O = _tau * _tau / 2;
        if (eq_data.approximation == "p1") {
            u[1][j] = eq_data.psi1(x) + _tau * eq_data.psi2(x) + eq_data.psi1_dir2(x) * O; // разложение по Тейлору
        } else if (eq_data.approximation == "p2") { // разложение по Тейлору с учётом уравнения
            u[1][j] = eq_data.psi1(x) + _tau * eq_data.psi2(x) +
                eq_data.a * O * eq_data.psi1_dir2(x) +
                eq_data.b * O * eq_data.psi1_dir1(x) +
                eq_data.c * O * eq_data.psi1(x) +
                eq_data.d * O * eq_data.psi2(x) +
                O * eq_data.f(x, 0);
        }
    }

    u[1][0] = eq_data.phi0(_tau);
    u[1][N - 1] = eq_data.phil(_tau);

    std::vector<double> a(N, 0);
    std::vector<double> b(N, 0);
    std::vector<double> c(N, 0);
    std::vector<double> d(N, 0);

    double damp_term = eq_data.d * _tau / 2.0;
    double conv_term = eq_data.b * _tau * _tau / (2.0 * _h);
    double reac_term = eq_data.c * _tau * _tau;

    for (int k = 2; k < K; ++k) {
        double t = k * _tau;
        for (int j = 1; j < N - 1; ++j) {
            double x = j * _h;
            a[j] = -_sigma;
            b[j] = (1.0 + damp_term + 2.0 * _sigma);
            c[j] = -_sigma;
            d[j] =
                (2.0 + reac_term) * u[k - 1][j] -
                (1.0 - damp_term) * u[k - 2][j] +
                conv_term * (u[k-1][j+1] - u[k-1][j-1]) +
                _tau * _tau * eq_data.f(x, (k-1)*_tau);
        }

        a[0] = 0.0;
        b[0] = 1.0;
        c[0] = 0.0;
        d[0] = b[0] * eq_data.phi0(k * _tau);

        a[N - 1] = 0.0;
        b[N - 1] = 1.0;
        c[N - 1] = 0.0;
        d[N - 1] = eq_data.phil(k * _tau);

        u[k] = tma(a, b, c, d);
    }

    for (int k = 0; k < K; ++k) {
        for (int j = 0; j < N; ++j) {
            fout << u[k][j] << " ";
        }
        fout << "\n";
    }
    fout << std::flush;
    fout.close();

    return "implicit";
}