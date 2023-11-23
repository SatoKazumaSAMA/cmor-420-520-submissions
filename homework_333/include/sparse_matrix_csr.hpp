#ifndef SPARSE_MATRIX_CSR_HPP
#define SPARSE_MATRIX_CSR_HPP

#include "abstract_matrix.hpp"

class SparseMatrixCSR : public AbstractMatrix {
public:
    SparseMatrixCSR(int m, int n, int nnz, const double *values, const int *col_indices, const int *row_indices);

    ~SparseMatrixCSR();

    const double *get_nzval() const { return _nzval; }
    const int *get_col_index() const { return _col_index; }
    const int *get_row_index() const { return _row_index; }

    double operator()(int row, int column) const override;

private:
    double *_nzval;
    int *_col_index;
    int *_row_index;
};

Vector operator*(const SparseMatrixCSR &A, const Vector &x);

#endif // SPARSE_MATRIX_CSR_HPP
