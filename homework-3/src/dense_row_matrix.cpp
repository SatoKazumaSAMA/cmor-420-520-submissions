#include "dense_row_matrix.hpp"
#include <stdexcept>

DenseRowMatrix::DenseRowMatrix(int m, int n) : AbstractMatrix(m, n) {
    _data = new double[m * n](); // Initialize with zeros
}

DenseRowMatrix::~DenseRowMatrix() {
    delete[] _data;
}

double DenseRowMatrix::operator()(int row, int column) const {
    if (row >= num_rows() || column >= num_columns()) {
        throw std::out_of_range("Index out of bounds");
    }
    return _data[row * num_columns() + column];
}

double &DenseRowMatrix::operator()(int row, int column) {
    if (row >= num_rows() || column >= num_columns()) {
        throw std::out_of_range("Index out of bounds");
    }
    return _data[row * num_columns() + column];
}
