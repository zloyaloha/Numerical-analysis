#include "tridiagonal-matrix.h"
#include <vector>

int main() {
    Tridiagonal tr;
    Matrix<double> mat("/home/zloyaloha/Yandex.Disk/MAI/Numerical-analysis/lab1/tridiagonal-matrix/test/test.txt");
    std::vector<double> res = tr.solve(mat);
    for (auto& x: res) {
        std::cout << x << ' ';
    }
    std::cout << std::endl;
    
    mat.checkEquality(res);

    // lu.decompose(mat, 1);
}