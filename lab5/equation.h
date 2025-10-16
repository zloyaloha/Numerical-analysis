#include "tridiagonal-matrix.h"
#include <cmath>
#include <fstream>
#include <functional>
#include <sstream>
#include <stdexcept>
enum BoundType {
    First,
    Second,
    Third
};

enum SolverType {
    Explicit,
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
private:
    EquationData eq_data;
    std::function<void(int,int)> solve_func;
    double _h, _tau, _sigma, _a;
    std::ofstream fout;
public:
    Parabolic(const EquationData& data, const SolverType& type);
    void solve(int N, double sigma, double T, double a);
private:
    void explicit_solve(int N, int K);
    void implicit_solve(int N, int K);
    void crank_nicolson_solve(int N, int K);
    void exact(int N, int K);
};

class ParabolicPostProcessor {
private:
    std::ifstream ifs;
public:
    ParabolicPostProcessor(const std::string& file_name);
    void graphics_u(double h, const std::string& name, int plot_steps);
};

inline Parabolic::Parabolic(const EquationData& data, const SolverType& type) : eq_data(data)
{
    switch (type)
    {
    case Explicit:
        solve_func = [this](int N, int K) { return this->explicit_solve(N, K); };
        fout.open("explicit.txt");
        break;
    case Implicit:
        solve_func = [this](int N, int K) { return this->implicit_solve(N, K); };
        fout.open("implicit.txt");
        std::cout << "Открыт файл: " << "implicit.txt" << std::endl;
        break;
    case Crank_Nicholson:
        solve_func = [this](int N, int K) { return this->crank_nicolson_solve(N, K); };
        fout.open("crank-nicolson.txt");
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

inline void Parabolic::solve(int N, double sigma, double T, double a)
{
    _sigma = sigma;
    _a = a;
    _h = eq_data.l / N;
    _tau = _sigma * _h * _h / a;
    int K = std::ceil(T / _tau);
    fout << K << ' ' << N + 1 << std::endl;
    solve_func(N, K);
}

inline void print_arr(const std::vector<double>& vec) {
    for (int j = 0; j < vec.size(); ++j) {
        std::cout << vec[j] << ' ';
    }
    std::cout << std::endl;

}

inline void Parabolic::explicit_solve(int N, int K)
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
            u_curr[N] = eq_data.phil(k*_tau);
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
}


inline void Parabolic::implicit_solve(int N, int K)
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

        if (eq_data.bound_type == First) {
            a[0] = 0; b[0] = -(1+2*_sigma); c[0] = _sigma;
            d[0] = -(u_prev[0] + _sigma*eq_data.phi0(k*_tau));

            a[N] = _sigma; b[N] = -(1+2*_sigma); c[N-1] = 0;
            d[N] = -(u_prev[N] + _sigma*eq_data.phil(k*_tau));
        } else if (eq_data.bound_type == Second) {
            a[0] = 0;
            b[0] = -(1 + 2 * _sigma);
            c[0] = _sigma;
            d[0] = -(u_prev[0] + _sigma * eq_data.phi0(k * _tau)) - _tau * eq_data.f(0, k * _tau);
            a[N] = _sigma;
            b[N] = -(1 + 2 * _sigma);
            c[N] = 0;
            d[N] = -(u_prev[N] + _sigma * eq_data.phil(k * _tau)) - _tau * eq_data.f((N) * _h, k * _tau);
        }

        u_curr = tma(a,b,c,d);
        u_prev = u_curr;

        for (int i = 0; i < u_curr.size(); ++i) {
            fout << u_curr[i] << ' ';
        }
        fout << '\n';
    }
    fout << std::endl;
    fout.close();
    return;
}


inline void Parabolic::crank_nicolson_solve(int N, int K)
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

        a[0] = 0; c[0] = 0; b[0] = 1;
        d[0] = eq_data.phi0(k * _tau);

        a[N] = 0; c[N] = 0; b[N] = 1;
        d[N] = eq_data.phil(k * _tau);

        u_curr = tma(a, b, c, d);

        for (int j = 0; j <= N; ++j)
            fout << u_curr[j] << ' ';
        fout << '\n';

        u_prev = u_curr;
    }

    fout << std::endl;
}

inline void Parabolic::exact(int N, int K)
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
}

inline ParabolicPostProcessor::ParabolicPostProcessor(const std::string& filename) : ifs(filename)
{
    std::cout << filename << std::endl;
    if (!ifs.is_open()) {
        throw std::runtime_error("can't open file");
    }
}

inline void ParabolicPostProcessor::graphics_u(double h,
                                               const std::string& name,
                                               int plot_steps)
{
    int K, N;
    ifs >> K >> N; // <-- проверь порядок с файлом!

    std::vector<std::vector<double>> u(K, std::vector<double>(N, 0.0));
    for (int i = 0; i < K; ++i)
        for (int j = 0; j < N; ++j)
            ifs >> u[i][j];

    // === читаем точное решение (если есть)
    std::ifstream analytic_ifs("exact.txt");
    std::vector<std::vector<double>> u_exact;
    bool has_exact = false;

    if (analytic_ifs) {
        has_exact = true;
        int K_e, N_e;
        analytic_ifs >> K_e >> N_e;
        u_exact.assign(K_e, std::vector<double>(N_e + 1, 0.0));
        for (int i = 0; i < K_e; ++i)
            for (int j = 0; j < N_e; ++j)
                analytic_ifs >> u_exact[i][j];
    }

    std::vector<int> time_indices = {
        std::max(0, int(K * 0.05)),
        std::max(0, int(K * 0.10)),
        std::max(0, int(K * 0.25))
    };

    // === создаём временные файлы для Gnuplot
    std::ofstream temp("temp_data.txt");
    for (int j = 0; j < N; ++j) {
        temp << j * h;
        for (int idx : time_indices)
            temp << " " << u[idx][j];
        temp << "\n";
    }
    temp.close();

    if (has_exact) {
        std::ofstream temp_exact("temp_exact.txt");
        for (int j = 0; j < N; ++j) {
            temp_exact << j * h;
            for (int idx : time_indices)
                temp_exact << " " << u_exact[idx][j];
            temp_exact << "\n";
        }
        temp_exact.close();
    }

    // === генерируем Gnuplot скрипт
    std::stringstream gnuplotCmd;
    gnuplotCmd << "set title '" << name << " over time'\n"
               << "set xlabel 'x'\n"
               << "set ylabel 'u(x,t)'\n"
               << "set grid\n"
               << "plot ";

    for (size_t i = 0; i < time_indices.size(); ++i) {
        gnuplotCmd << "'temp_data.txt' using 1:" << (i + 2)
                   << " with lines title 'num t=" << time_indices[i] << "'";
        if (i + 1 != time_indices.size() || has_exact)
            gnuplotCmd << ", \\\n     ";
    }

    if (has_exact) {
        for (size_t i = 0; i < time_indices.size(); ++i) {
            gnuplotCmd << "'temp_exact.txt' using 1:" << (i + 2)
                       << " with lines dashtype 2 title 'exact t=" << time_indices[i] << "'";
            if (i + 1 != time_indices.size())
                gnuplotCmd << ", \\\n     ";
        }
    }

    gnuplotCmd << "\npause -1\n";

    std::ofstream script("gnuplot_script.gp");
    script << gnuplotCmd.str();
    script.close();

    system("gnuplot gnuplot_script.gp");
}
