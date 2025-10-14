#include "qr.h"
#include <vector>

int main() {
    QR qr;
    double eps;
    std::cin >> eps;
    Matrix<double> mat("/home/zloyaloha/Yandex.Disk/MAI/Numerical-analysis/lab1/qr/test/test2.txt");
    qr.decompose(mat, eps);
}