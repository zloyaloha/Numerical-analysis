#include "tridiagonal-matrix.h"
#include <cmath>
#include <fstream>
#include <functional>
#include <stdexcept>
#include "DU_postprocessor/post_processor.h"

enum BoundType {
    First,
    Second,
    Third
};

enum SolverType {
    Explicit = 1,
    Implicit,
    Crank_Nicholson,
    Exact
};

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
    std::function<std::string(int,int)> solve_func;
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

inline Parabolic::Parabolic(const EquationData& data, const SolverType& type) : eq_data(data)
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

inline std::string Parabolic::solve(int N, double sigma, double T, double a)
{
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

inline std::string Parabolic::explicit_solve(int N, int K)
{
    std::vector<double> u_prev(N + 1, 0.0), u_curr(N + 1, 0.0);

    for (int j = 0; j < N; ++j) {
        u_prev[j] = eq_data.psi(j * _h);
    }
    for (int j = 0; j < u_prev.size(); ++j) {
        fout << u_prev[j] << ' ';
    }
    fout << '\n';

    for (int k = 1; k < K; ++k) {
        u_curr[0] = eq_data.phi0(k * _tau); // первый слой
        if (eq_data.bound_type == First) { // последний слой
            u_curr[N] = u_curr[N-2] + eq_data.phil(k*_tau) * _h;
        } else if (eq_data.bound_type == Second) {
            // u_curr[N] = eq_data.phil(k * _tau);
            u_curr[N] = u_curr[N-1] + _h * eq_data.phil(k * _tau);;
            // u_curr[N] = (4.0 * u_curr[N-1] - u_curr[N-2] + 2.0 * _h * eq_data.phil(k * _tau)) / 3.0;
        } else if (eq_data.bound_type == Third) {
            u_curr[N] = (eq_data.phil(k*_tau) + u_curr[N-2]/_h + 2*_tau*u_prev[N-1]/_h)
            / (1/_h + 2*_tau/_h);
        }

        for (int j = 1; j < N; ++j) {
            u_curr[j] = _sigma * u_prev[j+1]
                      + (1 - 2*_sigma) * u_prev[j]
                      + _sigma * u_prev[j-1]
                      + _tau * eq_data.f(j * _h, (k-1) * _tau);
            // std::cout << "k = " << k << " j= " << j << "_sigma = " << _sigma << " u_prev[j+1] = " << u_prev[j+1] << " u_prev[j] = " << u_prev[j] << " u_prev[j - 1] = " << u_prev[j - 1] << " u_curr[j]= " << u_curr[j] << std::endl;
        }

        for (int j = 0; j < u_curr.size(); ++j) {
            fout << u_curr[j] << ' ';
            // std::cout << u_curr.size() << std::endl;
        }
        fout << '\n';

        u_prev = u_curr;
    }
    fout << std::endl;
    fout.close();
    return "explicit";
}


inline std::string Parabolic::implicit_solve(int N, int K)
{
    std::vector<double> a(N+1,0), b(N+1,0), c(N+1,0), d(N+1, 0);
    std::vector<double> u_prev(N+1, 0.0), u_curr(N+1, 0.0);

    for (int i = 1; i < N; ++i)
        u_prev[i] = eq_data.psi(i*_h);

    for (int k = 1; k < K; ++k) {
        for (int j = 1; j < N; ++j) {
            a[j] = _sigma;
            b[j] = -(1 + 2*_sigma);
            c[j] = _sigma;
            d[j] = -u_prev[j] - _tau*eq_data.f(j*_h, k*_tau);
        }

        a[0] = 0;
        b[0] = -(1 + 2 * _sigma);
        c[0] = _sigma;
        d[0] = -(u_prev[0] + _sigma * eq_data.phi0(k * _tau)) - _tau * eq_data.f(0, k * _tau);
        a[N] = _sigma;
        b[N] = -(1 + _sigma);
        c[N] = 0;
        d[N] = -(u_prev[N] + _tau * eq_data.f((N)*_h, k*_tau)
                    + _sigma * _h * eq_data.phil(k*_tau));

        // a[0] = 0;
        // b[0] = 1.0;   // просто фиксируем u_0 = phi0
        // c[0] = 0;
        // d[0] = eq_data.phi0(k * _tau);

        // // Правая граница j=N (Neumann, второй род, 2-й порядок)
        // a[N] = -1.0;          // коэффициент при u_{N-2}
        // b[N] = 3.0;           // коэффициент при u_N
        // c[N] = -4.0;          // коэффициент при u_{N-1}
        // d[N] = 2.0 * _h * eq_data.phil(k * _tau);

        u_curr = tma(a,b,c,d);
        u_prev = u_curr;

        for (int i = 0; i < u_curr.size(); ++i) {
            fout << u_curr[i] << ' ';
        }
        fout << '\n';
    }
    fout << std::endl;
    fout.close();
    return "implicit";
}


inline std::string Parabolic::crank_nicolson_solve(int N, int K)
{
    std::vector<double> a(N + 1, 0), b(N + 1, 0), c(N + 1, 0), d(N + 1, 0);
    std::vector<double> u_prev(N + 1, 0.0), u_curr(N + 1, 0.0);

    for (int j = 0; j <= N; ++j)
        u_prev[j] = eq_data.psi(j * _h);

    for (int j = 0; j <= N; ++j)
        fout << u_prev[j] << ' ';
    fout << '\n';

    for (int k = 1; k < K; ++k) {
        double t_half = (k - 0.5) * _tau;

        for (int j = 1; j < N; ++j) {
            a[j] = -_sigma / 2.0;
            b[j] = 1.0 + _sigma;
            c[j] = -_sigma / 2.0;

            d[j] = (_sigma / 2.0) * u_prev[j - 1]
                  + (1.0 - _sigma) * u_prev[j]
                  + (_sigma / 2.0) * u_prev[j + 1]
                  + _tau * eq_data.f(j * _h, t_half);
        }

        // a[0] = 0; c[0] = 0; b[0] = 1;
        // d[0] = eq_data.phi0(k * _tau);
        a[0] = 0; b[0] = 1.0; c[0] = 0;
        d[0] = eq_data.phi0(k * _tau);

        // a[N] = 0; c[N] = 0; b[N] = 1;
        // // d[N] = eq_data.phil(k * _tau);
        // d[N] = u_curr[N-1] + _h * eq_data.phil(k * _tau);
        a[N] = -1.0;     // коэффициент при u_{N-2}
        b[N] = -3.0;     // коэффициент при u_N
        c[N] = 4.0;      // коэффициент при u_{N-1}
        d[N] = 2.0 * _h * eq_data.phil(k * _tau);

        u_curr = tma(a, b, c, d);

        for (int j = 0; j <= N; ++j)
            fout << u_curr[j] << ' ';
        fout << '\n';

        u_prev = u_curr;
    }

    fout << std::endl;
    return "crank-nicholson";
}

inline std::string Parabolic::exact(int N, int K)
{
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
