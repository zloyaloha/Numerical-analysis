#include <iostream>
#include <iterator>
#include <vector>
#include <cmath>
#include <algorithm>
#include "cubic-splines.h"

int main() {
    // std::vector<double> x = {0.1, 0.5, 0.9, 1.3, 1.7};
    // std::vector<double> y = {-2.2026, -0.19315, 0.79464, 1.5624, 2.2306};

    // std::vector<double> x = {-170, -130, -90, -50, -10};

    std::vector<double> x = {-1., 0., 3., 4.};
    std::vector<double> y = {-2., 6, 0, 1};

    double x_star = 0.8;
    CubicSplines spline(x, y);
    spline.validateSpline();

    std::cout << "S(" << x_star << ") = " << spline.evaluate(x_star) << std::endl;

    std::string dataFile = "spline_data.dat";
    spline.writeToFile(dataFile, 0.01);

    std::ofstream gp("plot.gnuplot");
    gp << "set title 'Cubic Spline Interpolation'\n";
    gp << "set grid\n";
    gp << "plot '" << dataFile << "' with lines title 'Spline', '-' with points pt 7 lc rgb 'red' title 'Data Points'\n";

    for (size_t i = 0; i < x.size(); ++i) {
        gp << x[i] << " " << y[i] << "\n";
    }
    gp << "e\n";
    gp.close();

    std::system("gnuplot -p plot.gnuplot");  // -p keeps window open

    spline.printPolynomials();

    return 0;
}