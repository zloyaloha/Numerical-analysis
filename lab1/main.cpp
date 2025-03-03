#include "lu.h"
#include "matrix.h"
#include <fstream>

int main() {
    LU lu;
    Matrix<double> mat("/home/zloyaloha/programming/Numerical-analysis/lab1/test/test.txt");
    std::vector<double> res = lu.solve(mat);
    for (auto& x: res) {
        std::cout << x << ' ';
    }
    std::cout << std::endl;
    lu.decompose(mat, 1);
}