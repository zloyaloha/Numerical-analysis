#include <iostream>
#include <iterator>
#include <vector>
#include <cmath>
#include <algorithm>
#include "cubic-splines.h"

int main()
{
    
    std::vector<double> x = {0.1, 0.5, 0.9, 1.3, 1.7};
    std::vector<double> y = {-2.2026, -0.19315, 0.79464, 1.5624, 2.2306};

    double x_star = 0.8;
    CubicSplines s(x, y);

    std::cout << s.evaluate(x_star) << std::endl;
}