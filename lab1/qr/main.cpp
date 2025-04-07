#include "qr.h"
#include <vector>

int main() {
    QR qr;
    Matrix<double> mat("/home/zloyaloha/programming/Numerical-analysis/lab1/qr/test/test.txt");
    auto A = qr.decompose(mat, 0.01);
    A.print();
    // std::vector<double> eigenvalues;
    // for (int i = 0; i < R.getRows(); ++i) {
    //     eigenvalues.push_back(R.at(i, i));
    // }

    // std::vector<Matrix<double>> eigenvectors;
    // for (int i = 0; i < Q.getCols(); ++i) {
    //     Matrix<double> eigenvector(Q.getCols(), 1, 0);
    //     for (int j = 0; j < Q.getRows(); ++j) {
    //         eigenvector.at(j, 0) = Q.at(j, i);
    //     }
    //     eigenvectors.push_back(eigenvector);
    //     std::cout << eigenvalues[i] << std::endl;
    //     eigenvectors[i].print();
    // }

    

    // for (int i = 0; i < eigenvalues.size(); ++i) {
    //     std::cout << "========Check " << i << "========" << std::endl;

    //     std::cout << "matrix * eigenvector" << std::endl;
    //     (mat * eigenvectors[i]).print();
    //     std::cout << std::endl;

    //     std::cout << "eigenvalue * eigenvector" << std::endl;
    //     (eigenvectors[i] * eigenvalues[i]).print();
    //     std::cout << std::endl;
    //     std::cout << "======================" << std::endl;
    // }
}