#pragma once
#include <iostream>
#include <fstream>
#include <vector>

template <typename T>
class Matrix {
private:
    std::vector<std::vector<T>> data;
    size_t rows;
    size_t cols;

public:
    Matrix(size_t rows, size_t cols) : rows(rows), cols(cols) {
        data.resize(rows, std::vector<T>(cols));
    }

    Matrix(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Ошибка открытия файла");
        }

        file >> rows >> cols;
        data.resize(rows, std::vector<T>(cols, 0.));

        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                file >> data[i][j];
            }
        }

        file.close();
    }

    Matrix(size_t rows_, size_t cols_, T value) : rows(rows_), cols(cols_), data(rows_, std::vector<T>(cols_, value)) {}

    T& at(size_t row, size_t col) {
        if (row >= rows || col >= cols) {
            throw std::out_of_range("Индекс " + std::to_string(row) + " " + std::to_string(col) +  " выходит за пределы матрицы");
        }
        return data[row][col];
    }

    const T& at(size_t row, size_t col) const {
        if (row >= rows || col >= cols) {
            throw std::out_of_range("Индекс выходит за пределы матрицы");
        }
        return data[row][col];
    }

    std::vector<T>& at(size_t row) {
        if (row >= rows) {
            throw std::out_of_range("Индекс выходит за пределы матрицы");
        }
        return data[row];
    }

    const std::vector<T>& at(size_t row) const {
        if (row >= rows) {
            throw std::out_of_range("Индекс выходит за пределы матрицы");
        }
        return data[row];
    }

    Matrix<T> multiply(const Matrix<T>& other) const {
        if (cols != other.rows) {
            throw std::invalid_argument("Число столбцов первой матрицы должно быть равно числу строк второй матрицы");
        }

        Matrix<T> result(rows, other.cols);

        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < other.cols; ++j) {
                result.at(i, j) = 0;
                for (size_t k = 0; k < cols; ++k) {
                    result.at(i, j) += data[i][k] * other.at(k, j);
                }
            }
        }

        return result;
    }

    void swapRows(size_t row1, size_t row2) {
        if (row1 >= rows || row2 >= rows) {
            throw std::out_of_range("Индекс строки выходит за пределы матрицы");
        }

        std::swap(data[row1], data[row2]);
    }

    void print() const {
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                std::cout << data[i][j] << "\t\t";
            }
            std::cout << std::endl;
        }
    }

    void checkEquality(const std::vector<double>& res)
    {
        double x_ = 0;
        for (int i = 0; i < this->getRows(); ++i) {
            for (int j = 0; j < this->getCols() - 1; ++j) {
                x_ += this->at(i, j) * res[j];
            }
            std::cout << x_ << ' ' << this->at(i, this->getCols() - 1) << std::endl;
            x_ = 0;
        }
    }

    size_t getRows() const { return rows; }
    size_t getCols() const { return cols; }
};