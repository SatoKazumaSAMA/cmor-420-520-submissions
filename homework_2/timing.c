#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "include/matrix.h" // Adjust the path as needed

double time_matvec(void (*matvec_function)(double *, const Matrix *, const double *), 
                   const Matrix *A, 
                   const double *x, 
                   double *b, 
                   int num_samples) {
    clock_t start, end;
    double cpu_time_used;

    start = clock();
    for (int i = 0; i < num_samples; ++i) {
        matvec_function(b, A, x);
    }
    end = clock();

    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    return cpu_time_used / num_samples;
}

int main(int argc, char *argv[]) {
    int m = 1000, n = 1000, num_samples = 10;

    if(argc == 4) {
        m = atoi(argv[1]);
        n = atoi(argv[2]);
        num_samples = atoi(argv[3]);
    }

    Matrix A_cont = init_matrix_contiguous(m, n);
    Matrix A_non_cont = init_matrix_non_contiguous(m, n);
    double x[n];
    double b[m];
    
    for (int j = 0; j < n; ++j) {
        x[j] = 1.0;
    }

    double avg_time = time_matvec(matvec, &A_cont, x, b, num_samples);
    printf("Average time for matvec with contiguous storage: %f seconds.\n", avg_time);

    avg_time = time_matvec(matvec, &A_non_cont, x, b, num_samples);
    printf("Average time for matvec with non-contiguous storage: %f seconds.\n", avg_time);

    avg_time = time_matvec(matvec_transpose, &A_cont, x, b, num_samples);
    printf("Average time for matvec_transpose with contiguous storage: %f seconds.\n", avg_time);

    avg_time = time_matvec(matvec_transpose, &A_non_cont, x, b, num_samples);
    printf("Average time for matvec_transpose with non-contiguous storage: %f seconds.\n", avg_time);

    free_matrix_contiguous(&A_cont);
    free_matrix_non_contiguous(&A_non_cont);
    
    return 0;
}
