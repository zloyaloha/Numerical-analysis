#include <vector>
#include <functional>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <sstream>

struct Point {
    double x, y1, y2;
};

class Solver {
public:
    virtual std::vector<Point> solve(std::function<double(const double&, const double&, const double&)> f, const double& h) = 0;
    void error(const std::vector<Point>& pt, const std::vector<Point>& pt_h2, std::function<double(const double&)> f, double coef);
    void grafics(const std::vector<Point>& pt, std::function<double(const double&)> f);
protected:
    std::string name;
    double p;
};

class CauchySolver : public Solver {
public:
    CauchySolver(double x0_, double y0_, double dy0_, double x_end_) : x0(x0_), y0(y0_), dy0(dy0_), x_end(x_end_) {}
protected:
    double x0, y0, dy0;
    double x_end;
};

class EulerSolver : public CauchySolver {
public:
    EulerSolver(double x0_, double y0_, double dy0_, double x_end_) : CauchySolver(x0_, y0_, dy0_, x_end_) {p = 1; name = "Euler";}
    std::vector<Point> solve(std::function<double(const double&, const double&, const double&)> f, const double& h) override;
};

class EulerCauchySolver : public CauchySolver {
public:
    EulerCauchySolver(double x0_, double y0_, double dy0_, double x_end_) : CauchySolver(x0_, y0_, dy0_, x_end_) {p = 2; name = "Euler-Cauchy";}
    std::vector<Point> solve(std::function<double(const double&, const double&, const double&)> f, const double& h);
};

class RungeKuttaSolver : public CauchySolver {
public:
    RungeKuttaSolver(double x0_, double y0_, double dy0_, double x_end_) : CauchySolver(x0_, y0_, dy0_, x_end_) {p = 4; name = "Runge-Kutta";}
    std::vector<Point> solve(std::function<double(const double&, const double&, const double&)> f, const double& h);
};

class AdamsSolver : public CauchySolver {
public:
    AdamsSolver(double x0_, double y0_, double dy0_, double x_end_) : CauchySolver(x0_, y0_, dy0_, x_end_) {p = 4; name = "Adams";}
    std::vector<Point> solve(std::function<double(const double&, const double&, const double&)> f, const double& h);
};

std::vector<Point> EulerSolver::solve(std::function<double(const double&, const double&, const double&)> f, const double& h)
{
    std::vector<Point> results;
    double x = x0;
    double y1 = y0;
    double y2 = dy0;

    while (x <= x_end + 1e-9) {
        results.push_back({x, y1, y2});

        double y1_new = y1 + h * y2;
        double y2_new = y2 + h * f(x, y1, y2);

        x += h;
        y1 = y1_new;
        y2 = y2_new;
    }

    return results;
}

std::vector<Point> EulerCauchySolver::solve(std::function<double(const double&, const double&, const double&)> f, const double& h)
{
    std::vector<Point> res;
    double x = x0;
    double y1 = y0;
    double y2 = dy0;

    while (x <= x_end + 1e-9) {
        res.push_back({x, y1, y2});

        double y1_pred = y1 + h * y2;
        double y2_pred = y2 + h * f(x, y1, y2);

        double f1_now = y2;
        double f2_now = f(x, y1, y2);
        double f1_pred = y2_pred;
        double f2_pred = f(x + h, y1_pred, y2_pred);

        y1 += (h / 2.0) * (f1_now + f1_pred);
        y2 += (h / 2.0) * (f2_now + f2_pred);
        x += h;
    }
    return res;
}

std::vector<Point> RungeKuttaSolver::solve(std::function<double(const double&, const double&, const double&)> f, const double& h)
{
    std::vector<Point> res;
    double x = x0;
    double y1 = y0;
    double y2 = dy0;

    while (x <= x_end + 1e-9) {
        res.push_back({x, y1, y2});

        double k0_1 = h * y2;
        double l0_1 = h * f(x, y1, y2);

        double k0_2 = h * (y2 + l0_1 / 2);
        double l0_2 = h * f(x + h / 2, y1 + k0_1 / 2, y2 + l0_1 / 2);

        double k0_3 = h * (y2 + l0_2 / 2);
        double l0_3 = h * f(x + h / 2, y1 + k0_2 / 2, y2 + l0_2 / 2);

        double k0_4 = h * (y2 + l0_3);
        double l0_4 = h * f(x + h, y1 + k0_3, y2 + l0_3);

        y1 += (1 / 6.0) * (k0_1 + 2 * k0_2 + 2 * k0_3 + k0_4);
        y2 += (1 / 6.0) * (l0_1 + 2 * l0_2 + 2 * l0_3 + l0_4);
        x += h;
    }
    return res;
}

std::vector<Point> AdamsSolver::solve(std::function<double(const double&, const double&, const double&)> f, const double& h)
{
    RungeKuttaSolver solver(x0, y0, dy0, x0 + h * 4);
    std::vector<Point> states = solver.solve(f, h);


    std::vector<Point> results(states);

    double x = states.back().x;
    while (x + h <= x_end + 1e-9) {
        Point& s0 = states[states.size() - 1];
        Point& s1 = states[states.size() - 2];
        Point& s2 = states[states.size() - 3];
        Point& s3 = states[states.size() - 4];

        double f1_0 = s0.y2;
        double f1_1 = s1.y2;
        double f1_2 = s2.y2;
        double f1_3 = s3.y2;

        double f2_0 = f(s0.x, s0.y1, s0.y2);
        double f2_1 = f(s1.x, s1.y1, s1.y2);
        double f2_2 = f(s2.x, s2.y1, s2.y2);
        double f2_3 = f(s3.x, s3.y1, s3.y2);

        double y1_new = s0.y1 + h / 24.0 * (55 * f1_0 - 59 * f1_1 + 37 * f1_2 - 9 * f1_3);
        double y2_new = s0.y2 + h / 24.0 * (55 * f2_0 - 59 * f2_1 + 37 * f2_2 - 9 * f2_3);
        x += h;

        results.push_back({x, y1_new, y2_new});
        states.push_back({x, y1_new, y2_new});
    }

    return results;
}

void Solver::grafics(const std::vector<Point>& results, std::function<double(const double&)> exact)
{

    std::ofstream temp("temp_data.txt");
    for (const auto& pt : results) {
        temp << pt.x << " " << pt.y1 << " " << exact(pt.x) << " " << pt.y2 << "\n";
    }
    temp.close();

    std::stringstream gnuplotCmd;
     gnuplotCmd << "set title '" << name << " vs Exact Solution'\n"
               << "set xlabel 'x'\n"
               << "set ylabel 'y'\n"
               << "set grid\n"
               << "plot 'temp_data.txt' using 1:2 with linespoints title '" << name << "', \\\n"
               << "     'temp_data.txt' using 1:3 with lines title 'Exact'\n"
               << "pause -1\n";

    std::ofstream script("gnuplot_script.gp");
    script << gnuplotCmd.str();
    script.close();

    system("gnuplot gnuplot_script.gp");
}

void Solver::error(const std::vector<Point>& results, const std::vector<Point>& results_h2, std::function<double(const double&)> exact, double coef)
{
    std::cout << '\n' << name << " table" << '\n';
    std::cout << std::fixed << std::setprecision(12);
    std::cout << "x\tY(h)\tY(h/2)\tRungeErr\tExact\tAbsErr" << std::endl;
    for (size_t i = 0; i < results.size(); ++i) {
        double x = results[i].x;
        double yh = results[i].y1;

        auto it = std::find_if(results_h2.begin(), results_h2.end(), [&](const Point& p) {
            return std::abs(p.x - x) < 1e-10;
        });

        if (it == results_h2.end()) continue;

        double yh2 = it->y1;

        double err_runge = std::fabs((yh2 - yh) / (std::pow(coef, p) - 1));

        double y_exact = exact(x);
        double abs_err = std::fabs(yh - y_exact);

        std::cout << x << "\t" << yh << "\t" << yh2 << "\t" << err_runge << "\t\t" << y_exact << "\t" << abs_err << std::endl;
    }
}