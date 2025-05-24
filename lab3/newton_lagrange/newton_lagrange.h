#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

class LagrangeInterpolator {
public:
    LagrangeInterpolator(const std::vector<double>& basePoints, const std::vector<double>& _Y) : X(basePoints), Y(_Y) {}

    double interpolate(const double& x) const {
        if (x > X[X.size() - 1] || x < X[0]) {
            throw std::range_error("X not in range");
        }
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

    void setNewPoints(const std::vector<double>& points, const std::vector<double>& Y_) {
        X = points;
        Y = Y_;
    }

    std::string getPolynomial() const {
        std::ostringstream oss;
        for (int i = 0; i < X.size(); ++i) {
            if (i != 0) oss << " + ";
            oss << std::setprecision(6) << Y[i];
            
            for (int j = 0; j < X.size(); ++j) {
                if (i != j) {
                    oss << " * (x - " << X[j] << ")/(" << X[i] << " - " << X[j] << ")";
                }
            }
        }
        return oss.str();
    }

    void printPolynomial() const {
        std::cout << "Lagrange interpolation polynomial:\n";
        std::cout << "P(x) = " << getPolynomial() << "\n";
    }

private:
    std::vector<double> X;
    std::vector<double> Y;
};

class NewtonInterpolator {
public:
    NewtonInterpolator(const std::vector<double>& basePoints, const std::vector<double>& Y_) : X(basePoints), coefs(Y_) {
        computeCoefficients();
    }

    void setNewPoints(const std::vector<double>& points, const std::vector<double>& Y_) {
        X = points;
        coefs = Y_;
        computeCoefficients();
    }

    double interpolate(const double& x) const {
        if (x > X[X.size() - 1] || x < X[0]) {
            std::cerr << "Invalid point" << std::endl;
            return -1.;
        }
        double result = coefs[0];
        double term = 1.0;
        for (int i = 1; i < X.size(); ++i) {
            term *= (x - X[i - 1]);
            result += coefs[i] * term;
        }
        return result;
    }

    void add_point(const double& new_x, const double& y_x) {
        double val = y_x;
        for (int i = X.size() - 1; i >= 0; --i) {
            val = (val - coefs[i]) / (new_x - X[i]);
        }
        X.push_back(new_x);
        coefs.push_back(val);
    }

    std::string getPolynomial() const {
        std::ostringstream oss;
        oss << std::setprecision(6) << coefs[0];
        
        for (int i = 1; i < X.size(); ++i) {
            oss << " + " << coefs[i];
            for (int j = 0; j < i; ++j) {
                oss << "*(x - " << X[j] << ")";
            }
        }
        return oss.str();
    }

    void printPolynomial() const {
        std::cout << "Newton interpolation polynomial:\n";
        std::cout << "P(x) = " << getPolynomial() << "\n";
    }

private:
    void computeCoefficients() {
        int n = X.size();
        for (int j = 1; j < n; ++j) {
            for (int i = n - 1; i >= j; --i) {
                coefs[i] = (coefs[i] - coefs[i - 1]) / (X[i] - X[i - j]);
            }
        }
    }

    std::vector<double> coefs;
    std::vector<double> X;
};