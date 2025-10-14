#include "rotating.h"
#include <vector>

int main() {
    Rotating rot(0.0000001);
    Matrix<long double> mat("/home/zloyaloha/Yandex.Disk/MAI/Numerical-analysis/lab1/rotating/test/text.txt");
    auto [eigenvalues, eigenvectors] =  rot.solve(mat);

    for (int i = 0; i < eigenvalues.size(); ++i) {
        Matrix<long double> x(mat.getRows(), 1, 0);
        for (int j = 0; j < x.getRows(); ++j) {
            x.at(j, 0) = eigenvectors.at(j, i);
        }


        std::cout << "========Check " << i << "========" << std::endl;

        std::cout << "matrix * eigenvector" << std::endl;
        (mat * x).print();
        std::cout << std::endl;

        std::cout << "eigenvalue * eigenvector" << std::endl;
        (x * eigenvalues[i]).print();
        std::cout << std::endl;
        std::cout << "======================" << std::endl;
    }
}