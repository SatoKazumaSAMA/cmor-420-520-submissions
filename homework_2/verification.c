#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "include/matrix.h" // Adjust the path as needed

// Function to compute the 1-norm of the error between two vectors
double compute_error(double *computed, double *exact, int size) {
    double error = 0.0;
    for (int i = 0; i < size; ++i) {
        error += fabs(computed[i] - exact[i]);
    }
    return error;
}

int main() {
    int m = 3, n = 4;

    // Initializing matrix with known values
    double A_values[3][4] = {
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {9, 10, 11, 12}
    };

    // Known vector x
    double x[4] = {1, 1, 1, 1};

    // Known results for b = A*x and b = transpose(A)*x
    double b_exact[3] = {10, 26, 42};
    double b_transpose_exact[4] = {15, 18, 21, 24};

    // Compute products for contiguous matrix storage
    Matrix A_cont = init_matrix_contiguous(m, n);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            A_cont.data[i][j] = A_values[i][j];
    
    double b_cont[3];
    matvec(b_cont, &A_cont, x);
    double error_cont = compute_error(b_cont, b_exact, m);
    printf("The error in computing b=Ax for matvec with contiguous matrix storage is %.5e\n", error_cont);

    double b_cont_transpose[4];
    matvec_transpose(b_cont_transpose, &A_cont, x);
    double error_cont_transpose = compute_error(b_cont_transpose, b_transpose_exact, n);
    printf("The error in computing b=A^Tx for matvec_transpose with contiguous matrix storage is %.5e\n", error_cont_transpose);

    free_matrix_contiguous(&A_cont);

    // Compute products for non-contiguous matrix storage
    Matrix A_non_cont = init_matrix_non_contiguous(m, n);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            A_non_cont.data[i][j] = A_values[i][j];

    double b_non_cont[3];
    matvec(b_non_cont, &A_non_cont, x);
    double error_non_cont = compute_error(b_non_cont, b_exact, m);
    printf("The error in computing b=Ax for matvec with non-contiguous matrix storage is %.5e\n", error_non_cont);

    double b_non_cont_transpose[4];
    matvec_transpose(b_non_cont_transpose, &A_non_cont, x);
    double error_non_cont_transpose = compute_error(b_non_cont_transpose, b_transpose_exact, n);
    printf("The error in computing b=A^Tx for matvec_transpose with non-contiguous matrix storage is %.5e\n", error_non_cont_transpose);

    free_matrix_non_contiguous(&A_non_cont);

    return 0;
}
