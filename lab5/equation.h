#include <cmath>
#include <fstream>
#include <functional>
#include <stdexcept>

#include "DU_postprocessor/post_processor.h"
#include "tridiagonal-matrix.h"

enum BoundType { First = 1, Second, Third };

enum SolverType { Explicit = 1, Implicit, Crank_Nicholson, Exact };

struct EquationData {
    std::function<double(double, double)> f;
    std::function<double(double)> psi;
    std::function<double(double)> phi0;
    std::function<double(double)> phil;
    std::function<double(double, double)> exact_solution;
    BoundType bound_type;
    double l;
};

class Parabolic {
public:
    EquationData eq_data;
    std::function<std::string(int, int)> solve_func;
    double _h, _tau, _sigma, _a;
    std::ofstream fout;

public:
    Parabolic(const EquationData& data, const SolverType& type);
    std::string solve(int N, double sigma, double T, double a);

private:
    std::string explicit_solve(int N, int K);
    std::string implicit_solve(int N, int K);
    std::string crank_nicolson_solve(int N, int K);
    std::string exact(int N, int K);
};

inline Parabolic::Parabolic(const EquationData& data, const SolverType& type) : eq_data(data) {
    switch (type) {
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
        case Crank_Nicholson:
            solve_func = [this](int N, int K) { return this->crank_nicolson_solve(N, K); };
            fout.open("result.txt");
            std::cout << "crank-nicholson" << std::endl;
            break;
        case Exact:
            solve_func = [this](int N, int K) { return this->exact(N, K); };
            fout.open("exact.txt");
            break;
        default:
            break;
    }

    if (!fout.is_open()) {
        throw std::runtime_error("can't open file");
    }
}

inline std::string Parabolic::solve(int N, double sigma, double T, double a) {
    _sigma = sigma;
    _a = a;
    _h = eq_data.l / N;
    _tau = _sigma * _h * _h / a;
    int K = std::ceil(T / _tau);
    fout << K << ' ' << N + 1 << std::endl;
    return solve_func(N, K);
}

inline void print_arr(const std::vector<double>& vec) {
    for (int j = 0; j < vec.size(); ++j) {
        std::cout << vec[j] << ' ';
    }
    std::cout << std::endl;
}

inline std::string Parabolic::explicit_solve(int N, int K) {
    std::vector<double> u_prev(N + 1, 0.0), u_curr(N + 1, 0.0);

    for (int j = 0; j <= N; ++j) {
        u_prev[j] = eq_data.psi(j * _h);
        fout << u_prev[j] << ' ';
    }
    fout << '\n';

    for (int k = 1; k < K; ++k) {
        double t = k * _tau;

        u_curr[0] = eq_data.phi0(t);

        for (int j = 1; j < N; ++j) {
            u_curr[j] = _sigma * u_prev[j + 1] + (1 - 2 * _sigma) * u_prev[j] + _sigma * u_prev[j - 1] +
                        _tau * eq_data.f(j * _h, (k - 1) * _tau);
        }

        double phi_l = eq_data.phil(t);

        switch (eq_data.bound_type) {
            case First:
                u_curr[N] = u_curr[N - 1] + _h * phi_l;
                break;

            case Second: {
                double f_val = eq_data.f(N * _h, (k - 1) * _tau);
                u_curr[N] = u_prev[N] + _tau * ((2.0 * _a / (_h * _h)) * (u_prev[N - 1] - u_prev[N] + _h * phi_l) + f_val);
                break;
            }

            case Third:
                u_curr[N] = (4.0 * u_curr[N - 1] - u_curr[N - 2] + 2.0 * _h * phi_l) / 3.0;
                break;
        }

        u_prev = u_curr;
        for (double val : u_curr) fout << val << ' ';
        fout << '\n';
    }

    fout.close();
    return "explicit";
}

inline std::string Parabolic::implicit_solve(int N, int K) {
    std::vector<double> a(N + 1), b(N + 1), c(N + 1), d(N + 1);
    std::vector<double> u_prev(N + 1), u_curr(N + 1);

    for (int i = 0; i <= N; ++i) u_prev[i] = eq_data.psi(i * _h);

    for (int k = 1; k < K; ++k) {
        double t = k * _tau;

        a[0] = 0;
        b[0] = 1;
        c[0] = 0;
        d[0] = eq_data.phi0(t);

        for (int j = 1; j < N; ++j) {
            a[j] = _sigma;
            b[j] = -(1 + 2 * _sigma);
            c[j] = _sigma;
            d[j] = -u_prev[j] - _tau * eq_data.f(j * _h, t);
        }

        double phi_l = eq_data.phil(t);

        switch (eq_data.bound_type) {
            case First:
                a[N] = -1.0;
                b[N] = 1.0;
                c[N] = 0.0;
                d[N] = _h * phi_l;
                break;

            case Second:
                a[N] = 2.0 * _sigma;
                b[N] = -(1.0 + 2.0 * _sigma);
                c[N] = 0.0;
                d[N] = -u_prev[N] - _tau * eq_data.f(N * _h, t) - 2.0 * _sigma * _h * phi_l;
                break;

            case Third:
                a[N] = 2.0 * _sigma;
                b[N] = -(1.0 + 2.0 * _sigma);
                c[N] = 0.0;
                d[N] = -u_prev[N] - _tau * eq_data.f(N * _h, t) - 2.0 * _sigma * _h * phi_l;

                break;
        }

        u_curr = tma(a, b, c, d);

        for (int i = 0; i < u_curr.size(); ++i) {
            fout << u_curr[i] << ' ';
        }
        fout << '\n';
        u_prev = u_curr;
    }
    fout << std::endl;
    fout.close();
    return "implicit";
}

inline std::string Parabolic::crank_nicolson_solve(int N, int K) {
    std::vector<double> a(N + 1), b(N + 1), c(N + 1), d(N + 1);
    std::vector<double> u_prev(N + 1, 0.0), u_curr(N + 1, 0.0);

    // Начальное условие
    for (int j = 0; j <= N; ++j) {
        u_prev[j] = eq_data.psi(j * _h);
        fout << u_prev[j] << ' ';
    }
    fout << '\n';

    for (int k = 1; k < K; ++k) {
        double t = k * _tau;  // Текущее время (k)
        double t_next = t;    // Время слоя k+1 (в цикле k - это слой, который мы ищем, значит t уже k*tau)
        double t_prev = (k - 1) * _tau;
        double t_half = (k - 0.5) * _tau;  // Время для функции источника f

        a[0] = 0;
        b[0] = 1;
        c[0] = 0;
        d[0] = eq_data.phi0(t);

        for (int j = 1; j < N; ++j) {
            a[j] = -_sigma / 2.0;
            b[j] = 1.0 + _sigma;
            c[j] = -_sigma / 2.0;

            d[j] = (_sigma / 2.0) * u_prev[j - 1] + (1.0 - _sigma) * u_prev[j] + (_sigma / 2.0) * u_prev[j + 1] +
                   _tau * eq_data.f(j * _h, t_half);
        }

        double phi_curr = eq_data.phil(t);       // Поток на слое k
        double phi_prev = eq_data.phil(t_prev);  // Поток на слое k-1

        switch (eq_data.bound_type) {
            case First:  // 1-й порядок: u_N - u_{N-1} = h * phi
                a[N] = -1.0;
                b[N] = 1.0;
                c[N] = 0.0;
                d[N] = _h * phi_curr;
                break;

            case Second:
            case Third:  // Используем Ghost Point для обоих случаев (для сохранения трехдиагональности)
            {
                a[N] = -_sigma;
                b[N] = 1.0 + _sigma;
                c[N] = 0.0;

                double f_val = eq_data.f(N * _h, t_half);

                d[N] = (1.0 - _sigma) * u_prev[N] + _sigma * u_prev[N - 1] + _sigma * _h * (phi_curr + phi_prev) + _tau * f_val;
                break;
            }
        }

        u_curr = tma(a, b, c, d);
        u_prev = u_curr;

        for (int j = 0; j <= N; ++j) fout << u_curr[j] << ' ';
        fout << '\n';
    }

    fout << std::endl;
    return "crank-nicholson";
}

inline std::string Parabolic::exact(int N, int K) {
    for (int k = 0; k < K; ++k) {
        for (int j = 0; j <= N; ++j) {
            double x = j * _h;
            double t = k * _tau;
            fout << eq_data.exact_solution(x, t) << " ";
        }
        fout << "\n";
    }

    fout << std::flush;
    return "exact";
}
