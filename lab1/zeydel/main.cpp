#include "zeydel.h"
#include <vector>

int main() {
    Zeydel tr;
    Matrix<double> mat("/home/zloyaloha/programming/Numerical-analysis/lab1/zeydel/test/test.txt");
    std::vector<double> res = tr.solve(mat);
    for (auto& x: res) {
        std::cout << x << ' ';
    }
    std::cout << std::endl;

    mat.checkEquality(res);

    // lu.decompose(mat, 1);
}