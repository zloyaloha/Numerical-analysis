#include "tridiagonal-matrix.h"
#include <functional>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>

std::vector<double> tma(const std::vector<double>& a,
                   const std::vector<double>& b,
                   const std::vector<double>& c,
                   const std::vector<double>& d)
{
    int n = b.size();
    std::vector<double> alpha(n), beta(n), x(n);

    alpha[0] = -c[0] / b[0];
    beta[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        double denom = b[i] + a[i] * alpha[i - 1];
        alpha[i] = -c[i] / denom;
        beta[i] = (d[i] - a[i] * beta[i - 1]) / denom;
    }

    x[n - 1] = beta[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }

    return x;
}

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
    std::function<std::vector<std::vector<double>>(int,int,double)> solve_func;
    double h, tau, sigma;
public:
    Parabolic(const EquationData& data, const SolverType& type);
    std::vector<std::vector<double>> solve(int N, int K, double T);
    void grafics_u(const std::vector<std::vector<double>>& u, double h, const std::string& name, int plot_steps = 5);
private:
    std::vector<std::vector<double>> explicit_solve(int N, int K, double T);
    std::vector<std::vector<double>> implicit_solve(int N, int K, double T);
    std::vector<std::vector<double>> crank_nicholson_solve(int N, int K, double T);
    std::vector<std::vector<double>> exact(int N, int K, double T);
};

Parabolic::Parabolic(const EquationData& data, const SolverType& type) : eq_data(data)
{
    switch (type)
    {
    case Explicit:
        solve_func = [this](int N, int K, double T) { return this->explicit_solve(N, K, T); };
        break;
    case Implicit:
        solve_func = [this](int N, int K, double T) { return this->implicit_solve(N, K, T); };
        break;
    case Crank_Nicholson:
        solve_func = [this](int N, int K, double T) { return this->crank_nicholson_solve(N, K, T); };
        break;
    case Exact:
        solve_func = [this](int N, int K, double T) { return this->exact(N, K, T); };
        break;
    default:
        break;
    }
}

std::vector<std::vector<double>> Parabolic::solve(int N, int K, double T)
{
    h = eq_data.l / N;
    tau = T / K;
    sigma = tau / h * h;
    return solve_func(N, K, T);
}

std::vector<std::vector<double>> Parabolic::explicit_solve(int N, int K, double T)
{
    std::vector<std::vector<double>> u(K, std::vector<double>(N, 0.0)); // сетка

    for (int j = 1; j < N - 1; ++j) {
        u[0][j] = eq_data.psi(j * h); // это мы знаем по условию?
    }

    for (int k = 1; k < K; ++k) {
        u[k][0] = eq_data.phi0(k * tau); // нижний слой

        for (int j = 1; j < N-1; ++j) {
            u[k][j] = sigma * u[k-1][j+1]
                    + (1 - 2 * sigma) * u[k-1][j]
                    + sigma * u[k-1][j-1]
                    + tau * eq_data.f(j*h, k*tau); // явная схема
        }

        if (eq_data.bound_type == First) { // граничные условия
            u[k][N-1] = u[k][N-2] + eq_data.phil(k*tau) * h;
        } else if (eq_data.bound_type == Second) {
            u[k][N-1] = eq_data.phil(k*tau);
        } else if (eq_data.bound_type == Third) {
            u[k][N-1] = (eq_data.phil(k*tau) + u[k][N-2]/h + 2*tau*u[k-1][N-1]/h)
                        / (1/h + 2*tau/h);
        }
    }
    return u;
}

std::vector<std::vector<double>> Parabolic::implicit_solve(int N, int K, double T)
{
    std::vector<double> a(N,0), b(N,0), c(N,0), d(N,0);
    std::vector<std::vector<double>> u(K, std::vector<double>(N, 0.0));

    for (int i = 1; i < N-1; ++i)
        u[0][i] = eq_data.psi(i*h);

    for (int k = 1; k < K; ++k) {
        for (int j = 1; j < N-1; ++j) {
            a[j] = sigma;
            b[j] = -(1 + 2*sigma);
            c[j] = sigma;
            d[j] = -u[k-1][j] - tau*eq_data.f(j*h, k*tau);
        }

        // примеры для граничных условий (a1p1)
        a[0] = 0; b[0] = -(1+2*sigma); c[0] = sigma;
        d[0] = -(u[k-1][0] + sigma*eq_data.phi0(k*tau));

        a[N-1] = sigma; b[N-1] = -(1+2*sigma); c[N-1] = 0;
        d[N-1] = -(u[k-1][N-1] + sigma*eq_data.phil(k*tau));

        u[k] = tma(a,b,c,d);
    }
    return u;
}


std::vector<std::vector<double>> Parabolic::crank_nicholson_solve(int N, int K, double T)
{
}

std::vector<std::vector<double>> Parabolic::exact(int N, int K, double T)
{
    std::vector<std::vector<double>> u(K, std::vector<double>(N, 0.0));
    for (int j = 0; j < N; ++j) {
        for (int k = 0; k < K; ++k) {
            u[k][j] = eq_data.exact_solution(j * h, k * tau);
        }
    }
    return u;
}

void Parabolic::grafics_u(const std::vector<std::vector<double>>& u, double h, const std::string& name, int plot_steps) {
    int K = u.size();
    if (K == 0) return;
    int N = u[0].size();

    std::ofstream temp("temp_data.txt");
    std::vector<int> time_indices;
    for (int i = 0; i < plot_steps; ++i) {
        int idx = i * (K - 1) / (plot_steps - 1);
        time_indices.push_back(idx);
    }

    for (int j = 0; j < N; ++j) {
        temp << j*h; // x
        for (int idx : time_indices) {
            temp << " " << u[idx][j]; // u в разные моменты времени
        }
        temp << "\n";
    }
    temp.close();

    // 2. Создаём скрипт gnuplot
    std::stringstream gnuplotCmd;
    gnuplotCmd << "set title '" << name << " over time'\n"
                << "set xlabel 'x'\n"
                << "set ylabel 'u'\n"
                << "set grid\n"
                << "plot ";

    for (size_t i = 0; i < time_indices.size(); ++i) {
        gnuplotCmd << "'temp_data.txt' using 1:" << (i + 2)
                    << " with lines title 't=" << time_indices[i] << "'";
        if (i + 1 != time_indices.size()) gnuplotCmd << ", \\\n     ";
    }
    gnuplotCmd << "\npause -1\n";

    std::ofstream script("gnuplot_script.gp");
    script << gnuplotCmd.str();
    script.close();

    // 3. Запускаем gnuplot
    system("gnuplot gnuplot_script.gp");
}