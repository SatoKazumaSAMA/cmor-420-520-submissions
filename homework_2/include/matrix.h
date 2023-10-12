#ifndef MATRIX_H
#define MATRIX_H

typedef struct {
    int m; // Number of rows
    int n; // Number of columns
    double **data; // Pointer to the 2D array of matrix entries
} Matrix;

Matrix init_matrix_contiguous(int m, int n);
Matrix init_matrix_non_contiguous(int m, int n);
void free_matrix_contiguous(Matrix *mat);
void free_matrix_non_contiguous(Matrix *mat);
void matvec(double *b, const Matrix *A, const double *x);
void matvec_transpose(double *b, const Matrix *A, const double *x);

#endif // MATRIX_H
