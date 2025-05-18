#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>

class LagrangeInterpolator {
public:
    LagrangeInterpolator(const std::vector<double>& basePoints, const std::vector<double>& _Y) : X(basePoints), Y(_Y) {}

    double interpolate(const double& x) const
    {
        double res = 0;
        for (int i = 0; i < X.size(); ++i) {
            double term = Y[i];
            for (int j = 0; j < X.size(); ++j) {
                if (i != j) {
                    term *= (x - X[j]) / (X[i] - X[j]);
                }
            }
            res += term;
        }
        return res;
    }

    void setNewPoints(const std::vector<double>& points, const std::vector<double>& Y_)
    {
        X = points;
        Y = Y_;
    }

private:
    std::vector<double> X;
    std::vector<double> Y;
};

class NewtonInterpolator {
public:
    NewtonInterpolator(const std::vector<double>& basePoints, const std::vector<double>& Y_) : X(basePoints), coefs(Y_)
    {
        int n = X.size();
        for (int j = 1; j < n; ++j) {
            for (int i = n - 1; i >= j; --i) {
                coefs[i] = (coefs[i] - coefs[i - 1]) / (X[i] - X[i - j]);
            }
        }
    }

    void setNewPoints(const std::vector<double>& points, const std::vector<double>& Y_)
    {
        X = points;
        coefs = Y_;
        int n = X.size();
        for (int j = 1; j < n; ++j) {
            for (int i = n - 1; i >= j; --i) {
                coefs[i] = (coefs[i] - coefs[i - 1]) / (X[i] - X[i - j]);
            }
        }
    }

    double interpolate(const double& x) const
    {
        double result = coefs[0];
        double term = 1.0;
        for (int i = 1; i < X.size(); ++i) {
            term *= (x - X[i - 1]);
            result += coefs[i] * term;
        }
        return result;
    }

    void add_point(const double& new_x, const double& y_x)
    {
        double val = y_x;
        for (int i = X.size() - 1; i >= 0; --i) {
            val = (val - coefs[i]) / (new_x - X[i]);
        }
        X.push_back(new_x);
        coefs.push_back(val);
    }

private:
    std::vector<double> coefs;
    std::vector<double> X;
};