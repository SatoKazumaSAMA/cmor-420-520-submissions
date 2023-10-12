#include <stdlib.h>
#include "matrix.h" // Include the header file to have access to the Matrix struct and function prototypes

// Function to initialize a matrix with contiguous memory layout
Matrix init_matrix_contiguous(int m, int n) {
    Matrix matrix;
    matrix.m = m;
    matrix.n = n;
    matrix.data = (double **)malloc(m * sizeof(double *));
    matrix.data[0] = (double *)malloc(m * n * sizeof(double));
    for(int i = 1; i < m; i++) {
        matrix.data[i] = matrix.data[i - 1] + n;
    }
    return matrix;
}

// Function to initialize a matrix with non-contiguous memory layout
Matrix init_matrix_non_contiguous(int m, int n) {
    Matrix matrix;
    matrix.m = m;
    matrix.n = n;
    matrix.data = (double **)malloc(m * sizeof(double *));
    for(int i = 0; i < m; i++) {
        matrix.data[i] = (double *)malloc(n * sizeof(double));
    }
    return matrix;
}

// Function to free the allocated memory for a matrix with contiguous memory layout
void free_matrix_contiguous(Matrix *mat) {
    free(mat->data[0]);
    free(mat->data);
}

// Function to free the allocated memory for a matrix with non-contiguous memory layout
void free_matrix_non_contiguous(Matrix *mat) {
    for(int i = 0; i < mat->m; i++) {
        free(mat->data[i]);
    }
    free(mat->data);
}

// Function to compute the matrix-vector product b = A*x
void matvec(double *b, const Matrix *A, const double *x) {
    for(int i = 0; i < A->m; i++) {
        b[i] = 0;
        for(int j = 0; j < A->n; j++) {
            b[i] += A->data[i][j] * x[j];
        }
    }
}

// Function to compute the transpose matrix-vector product b = transpose(A)*x
void matvec_transpose(double *b, const Matrix *A, const double *x) {
    for(int i = 0; i < A->n; i++) {
        b[i] = 0;
        for(int j = 0; j < A->m; j++) {
            b[i] += A->data[j][i] * x[j];
        }
    }
}

