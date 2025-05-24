#include <iostream>
#include <iterator>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include "newton_lagrange.h"

double f(const double& x)
{
    return std::log(x) + x;
}

void writeGraphs(const std::vector<double>& X,
                 const LagrangeInterpolator& lagrange,
                 const NewtonInterpolator& newton,
                 const std::vector<double>& Y,
                 double xmin, double xmax)
{
    std::ofstream dataFile("data.txt");
    for (size_t i = 0; i < X.size(); ++i) {
        dataFile << X[i] << " " << Y[i] << "\n";
    }
    dataFile.close();

    std::ofstream fitFile("fit.txt");
    for (double x = xmin; x <= xmax; x += 0.01) {
        double yL = lagrange.interpolate(x);
        double yN = newton.interpolate(x);
        fitFile << x << " " << yL << " " << yN << "\n";
    }
    fitFile.close();

    std::ofstream gp("plot.gp");
    gp << "set title 'Интерполяция'\n";
    gp << "set xlabel 'x'\n";
    gp << "set ylabel 'y'\n";
    gp << "plot \\\n";
    gp << "'data.txt' with points pt 7 lc rgb 'black' title 'Данные', \\\n";
    gp << "'fit.txt' using 1:2 with lines lw 2 lc rgb 'blue' title 'Лагранж', \\\n";
    gp << "'fit.txt' using 1:3 with lines lw 2 lc rgb 'red' title 'Ньютон'\n";
    gp << "pause -1\n";
    gp.close();

    system("gnuplot -persist plot.gp");
}

int main()
{
    std::vector<double> X2 = {0.1, 0.5, 0.9, 1.3};
    std::vector<double> Y2(X2.size());
    std::transform(X2.begin(), X2.end(), Y2.begin(), f);

    // std::vector<double> X1 = {0.1, 0.5, 1.1, 1.3};
    // std::vector<double> Y1(X1.size());
    // std::transform(X1.begin(), X1.end(), Y1.begin(), f);

    // std::vector<double> X1 = {10, 50, 90, 130};
    // std::vector<double> Y1(X1.size());
    // std::transform(X1.begin(), X1.end(), Y1.begin(), f);

    std::vector<double> X1 = {0., 3., 4.};
    std::vector<double> Y1 = {6, 0, 1};

    double x = -10;
    double real_val = f(x);

    LagrangeInterpolator inter (X1, Y1);
    NewtonInterpolator newton_inter(X1, Y1);

    double val = inter.interpolate(x);

    std::cout << "Приближение в точке x = " << x << " методом Лагранжа P(x) = " << val << '\n';
    std::cout  << "Реальное значение функции в точке f(x) = " << real_val << std::endl;
    std::cout << "Погрешность " << std::fabs(val - real_val) << std::endl;

    val = newton_inter.interpolate(x);
    newton_inter.printPolynomial();

    std::cout << "Приближение в точке x = " << x << " методом Ньютона P(x) = " << val << '\n';
    std::cout  << "Реальное значение функции в точке f(x) = " << real_val << std::endl;
    std::cout << "Погрешность " << std::fabs(val - real_val) << std::endl;

    double xmin = *std::min_element(X1.begin(), X1.end());
    double xmax = *std::max_element(X1.begin(), X1.end());

    writeGraphs(X1, inter, newton_inter, Y1, xmin, xmax);

    inter.setNewPoints(X2, Y2);
    val = inter.interpolate(x);
    inter.printPolynomial();

    std::cout << "Приближение в точке x = " << x << " методом Лагранжа P(x) = " << val << '\n';
    std::cout  << "Реальное значение функции в точке f(x) = " << real_val << std::endl;
    std::cout << "Погрешность " << std::fabs(val - real_val) << std::endl;

    newton_inter.setNewPoints(X2, Y2);
    val = newton_inter.interpolate(x);

    std::cout << "Приближение в точке x = " << x << " методом Ньютона P(x) = " << val << '\n';
    std::cout  << "Реальное значение функции в точке f(x) = " << real_val << std::endl;
    std::cout << "Погрешность " << std::fabs(val - real_val) << std::endl;

    xmin = *std::min_element(X2.begin(), X2.end());
    xmax = *std::max_element(X2.begin(), X2.end());

    writeGraphs(X2, inter, newton_inter, Y2, xmin, xmax);
}