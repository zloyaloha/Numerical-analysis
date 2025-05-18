#include <iostream>
#include <iterator>
#include <vector>
#include <cmath>
#include <algorithm>
#include "mnk.h"

int main() {
    MNK mnk;
    // std::vector<double> x = {0.1, 0.5, 0.9, 1.3, 1.7, 2.1};
    std::vector<double> y = {-2.2026, -0.19315, 0.79464, 1.5624, 2.2306, 2.8419};

    std::vector<double> x = {-210, -170, -130, -90, -50, -10};

    // std::vector<double> x = {-1., 0., 3., 4.};
    // std::vector<double> y = {-2., 6, 0, 1};

    std::vector<double> p1 = mnk.findCoefs(x, y, 1);
    std::vector<double> p2 = mnk.findCoefs(x, y, 2);
    std::vector<double> p3 = mnk.findCoefs(x, y, 3);

    std::cout << "Полином 1-й степени: y = " << p1[0] << " + " << p1[1] << "x\n";
    std::cout << "Ошибка " << mnk.calculateSSE(x, y, p1) << std::endl;
    std::cout << "Полином 2-й степени: y = " << p2[0] << " + " << p2[1] << "x + " << p2[2] << "x^2\n";
    std::cout << "Ошибка " << mnk.calculateSSE(x, y, p2) << std::endl;
    std::cout << "Полином 3-й степени: y = " << p3[0] << " + " << p3[1] << "x + " << p3[2] << "x^2\n";
    std::cout << "Ошибка " << mnk.calculateSSE(x, y, p3) << std::endl;

    std::ofstream data("data.txt");
    for (size_t i = 0; i < x.size(); ++i)
        data << x[i] << " " << y[i] << "\n";
    data.close();

    std::ofstream fit("fit.txt");
    double min = *std::min_element(x.begin(), x.end());
    double max = *std::max_element(x.begin(), x.end());
    for (double xi = min; xi <= max; xi += 0.1)
        fit << xi << " " << mnk.evaluatePoly(p1, xi) << " " << mnk.evaluatePoly(p2, xi) << " " << mnk.evaluatePoly(p3, xi) << "\n";
    fit.close();

    std::ofstream gp("plot.gp");
    gp << "set title 'Аппроксимация МНК'\n";
    gp << "set xlabel 'x'\n";
    gp << "set ylabel 'y'\n";
    gp << "plot \\\n";
    gp << "'data.txt' with points pt 7 lc rgb 'black' title 'Данные', \\\n";
    gp << "'fit.txt' using 1:2 with lines lw 2 lc rgb 'blue' title '1-я степень', \\\n";
    gp << "'fit.txt' using 1:3 with lines lw 2 lc rgb 'red' title '2-я степень', \\\n";
    gp << "'fit.txt' using 1:4 with lines lw 2 lc rgb 'green' title '3-я степень'\n";
    gp << "pause -1\n";
    gp.close();

    system("gnuplot -persist plot.gp");
    return 0;
}