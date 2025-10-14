#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>
#include "lu.h"
#include <cmath>

class MNK {
public:
    std::vector<double> findCoefs(const std::vector<double>& x, const std::vector<double>& y, int degree)
    {
        int n = x.size();
        int m = degree + 1;
        Matrix<double> A(m, m + 1, 0.);

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; ++j) {
                for (int k = 0; k < n; ++k)
                    A.at(i, j) += pow(x[k], i + j);
            }
            for (int k = 0; k < n; ++k)
                A.at(i, m) += y[k] * pow(x[k], i);
        }
        LU lu;
        return lu.solve(A);  
    }

    double evaluatePoly(const std::vector<double>& coeffs, double x)
    {
        double res = 0.0;
        for (int i = 0; i < coeffs.size(); ++i)
            res += coeffs[i] * pow(x, i);
        return res;
    }
    
    double calculateSSE(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& coeffs)
    {
        double sse = 0.0;
        for (int i = 0; i < x.size(); ++i) {
            double diff = y[i] - evaluatePoly(coeffs, x[i]);
            sse += diff * diff;
        }
        return sse;
    }

};