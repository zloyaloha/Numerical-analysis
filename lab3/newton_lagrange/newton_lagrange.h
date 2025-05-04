#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>

class I_Interpolator {
public:
    virtual ~I_Interpolator() = default;
    I_Interpolator(std::function<double(const double&)> f) : func(f) {}
    virtual double interpolate(const double& x) const = 0;
protected:
    std::function<double(const double&)> func;
};

class LagrangeInterpolator : public I_Interpolator {
public:
    LagrangeInterpolator(std::function<double(const double&)> f, const std::vector<double>& basePoints) : I_Interpolator(f), X(basePoints) {}

    double interpolate(const double& x) const override
    {
        std::vector<double> Y(X.size());
        std::transform(X.begin(), X.end(), Y.begin(), func);
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

    void setNewPoints(const std::vector<double>& points)
    {
        X = points;
    }

private:
    std::vector<double> X;
};

class NewtonInterpolator : public I_Interpolator {
public:
    NewtonInterpolator(std::function<double(const double&)> f, const std::vector<double>& basePoints) : I_Interpolator(f), X(basePoints), coefs(basePoints.size())
    {
        std::transform(X.begin(), X.end(), coefs.begin(), func);
        for (int j = 1; j < X.size(); ++j) {
            for (int i = X.size() - 1; i >= j; --i) {
                coefs[i] = (coefs[i] - coefs[i - 1]) / (X[i] - X[i - j]);
            }
        }
    }

    double interpolate(const double& x) const override
    {
        double result = coefs[0];
        double term = 1.0;
        for (int i = 1; i < X.size(); ++i) {
            term *= (x - X[i - 1]);
            result += coefs[i] * term;
        }
        return result;
    }

    void add_point(const double& new_x)
    {
        X.push_back(new_x);
        double y = func(new_x);
        for (int j = 1; j < X.size(); ++j) {
            y = (y - coefs[X.size() - j - 1]) / (new_x - X[X.size() - j - 1]);
        }
        coefs.push_back(y);
    }

    void setNewPoints(const std::vector<double>& points)
    {
        std::transform(X.begin(), X.end(), coefs.begin(), func);
        for (int j = 1; j < X.size(); ++j) {
            for (int i = X.size() - 1; i >= j; --i) {
                coefs[i] = (coefs[i] - coefs[i - 1]) / (X[i] - X[i - j]);
            }
        }
    }

private:
    std::vector<double> coefs;
    std::vector<double> X;
};