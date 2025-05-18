#include <iostream>
#include <iterator>
#include <vector>
#include <cmath>
#include <algorithm>
#include "integral.h"

double f(double x)
{
    return (std::pow(x, 2) / (std::pow(x, 4 ) + 256));
}

int main() {
    Integral integral(0., 2.5);

    double I_rect_1 = integral.findIntegralRectangle(f, 0.25);
    double I_rect_2 = integral.findIntegralRectangle(f, 0.5);
    double I_rect_RRR = integral.rungeRomberg(I_rect_2, I_rect_1, 2);

    double I_trapezoid_1 = integral.findIntegralTrapezoid(f, 0.25);
    double I_trapezoid_2 = integral.findIntegralTrapezoid(f, 0.5);
    double I_trapezoid_RRR = integral.rungeRomberg(I_trapezoid_2, I_trapezoid_1, 2);

    double I_simpson_1 = integral.findIntegralSimpson(f, 0.25);
    double I_simpson_2 = integral.findIntegralSimpson(f, 0.5);
    double I_simpson_RRR = integral.rungeRomberg(I_simpson_2, I_simpson_1, 4);


    std::cout << "Прямоугольники с шагом 0.25: " << I_rect_1 << std::endl;
    std::cout << "Прямоугольники с шагом 0.5: " << I_rect_2 << std::endl;
    std::cout << "Прямоугольники уточненные: " << I_rect_RRR << std::endl;

    std::cout << "Трапеция с шагом 0.25: " << I_trapezoid_1 << std::endl;
    std::cout << "Трапеция с шагом 0.5: " << I_trapezoid_2 << std::endl;
    std::cout << "Трапеция уточненная: " << I_trapezoid_RRR << std::endl;

    std::cout << "Симпсон с шагом 0.25: " << I_simpson_1 << std::endl;
    std::cout << "Симпсон с шагом 0.5: " << I_simpson_2 << std::endl;
    std::cout << "Симпсон уточненный: " << I_simpson_RRR << std::endl;
}