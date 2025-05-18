#include <iostream>
#include <iterator>
#include <vector>
#include <cmath>
#include <algorithm>
#include "derivative.h"

int main() {
    Derivative der;
    // std::vector<double> x = {0., 1, 2, 3, 4};
    // std::vector<double> y = {0, 2.0, 3.4142, 4.7321, 6.0};
    std::vector<double> x = {0., 0.1, 0.2, 0.3, 0.4};
    std::vector<double> y = {1., 1.1052, 1.2214, 1.3499, 1.4918};

   std::cout << der.findFirstDerivative(x, y, 0.25) << std::endl;
   std::cout << der.findSecondDerivative(x, y, 0.25) << std::endl;
}