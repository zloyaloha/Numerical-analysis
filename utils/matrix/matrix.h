#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <functional>

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

    Matrix(const std::vector<T>& main, const std::vector<T>& down, const std::vector<T>& upper, const std::vector<T>& right) : rows(main.size()), cols(main.size() + 1) {
        int n = main.size();
        data.resize(n, std::vector<T>(n + 1, 0.));

        for (int i = 0; i < n; ++i) {
            data[i][i] = main[i];

            if (i > 0) {
                data[i][i - 1] = down[i - 1];
            }

            if (i < n - 1) {
                data[i][i + 1] = upper[i];
            }

            data[i][n] = right[i];
        }
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

    Matrix<T> operator*(const Matrix<T>& other) const {
        if (cols != other.rows) {
            throw std::invalid_argument("Число столбцов первой матрицы должно быть равно числу строк второй матрицы");
        }
        Matrix<T> result(rows, other.cols, 0);  // Результирующая матрица, инициализирована нулями
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < other.cols; ++j) {
                for (size_t k = 0; k < cols; ++k) {
                    result.data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }
        return result;
    }

    Matrix<T> operator-(const Matrix<T>& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Размеры матриц должны быть одинаковыми для вычитания.");
        }
        Matrix<T> result(rows, cols, 0);  // Обратите внимание на использование cols, а не other.cols
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {  // Используем cols, так как оба массива должны быть одинаковыми по размеру
                result.data[i][j] = data[i][j] - other.data[i][j];
            }
        }

        return result;
    }

    Matrix<T> operator*(const T& alpha) const {

        Matrix<T> result(rows, cols);

        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result.at(i, j) = data[i][j] * alpha;
            }
        }

        return result;
    }

    Matrix<T> operator~() const {
        Matrix<T> result(cols, rows);  // Создаем новую матрицу с обратными размерами
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result.data[j][i] = data[i][j];  // Копируем элементы с транспонированием
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

    static Matrix<T> identity(int n) {
        Matrix<T> I(n, n, 0);  // Инициализируем матрицу, заполняя её нулями
        for (int i = 0; i < n; ++i) {
            I.data[i][i] = static_cast<T>(1);  // Присваиваем единицу на диагонали
        }
        return I;
    }

    size_t getRows() const { return rows; }
    size_t getCols() const { return cols; }
};