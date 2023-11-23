#include "sparse_matrix_csr.hpp"
SparseMatrixCSR::SparseMatrixCSR(int m, int n, int nnz, const double *values, const int *col_indices, const int *row_indices)
    : AbstractMatrix(m, n), _nzval(nullptr), _col_index(nullptr), _row_index(nullptr) {
    _nzval = new double[nnz];
    _row_index = new int[m + 1];
    _col_index = new int[nnz];

    for (int i = 0; i < nnz; ++i) {
        _nzval[i] = values[i];
        _col_index[i] = col_indices[i];
    }

    for (int i = 0; i <= m; ++i) {
        _row_index[i] = row_indices[i];
    }
}

SparseMatrixCSR::~SparseMatrixCSR()
{
    delete[] _nzval;
    delete[] _col_index;
    delete[] _row_index;
}


double SparseMatrixCSR::operator()(int row, int column) const {
    int row_start = _row_index[row];
    int row_end = _row_index[row + 1];
    for (int i = row_start; i < row_end; ++i) {
        if (_col_index[i] == column) {
            return _nzval[i];
        }
    }
    return 0.0; // Default value for missing entries
}


Vector operator*(const SparseMatrixCSR& A, const Vector& x) {
    const double* nzval = A.get_nzval();
    const int* row_index = A.get_row_index();
    const int* col_index = A.get_col_index();

    Vector result(A.num_rows());
    for (int i = 0; i < A.num_rows(); ++i) {
        int row_start = row_index[i];
        int row_end = row_index[i + 1];
        double sum = 0.0;
        for (int j = row_start; j < row_end; ++j) {
            int col = col_index[j];
            sum += nzval[j] * x(col);
        }
        result(i) = sum;
    }
    return result;
}
