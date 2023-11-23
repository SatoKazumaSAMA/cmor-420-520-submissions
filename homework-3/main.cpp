#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include "dense_row_matrix.hpp"
#include "sparse_matrix_csr.hpp"
#include "vector.hpp"

// Define the norm function to calculate the Euclidean norm of a Vector
double norm(const Vector& v) {
    double sum = 0.0;
    for (int i = 0; i < v.length(); ++i) {
        sum += v(i) * v(i);
    }
    return std::sqrt(sum);
}

int main() {
    const int n = 100;
    const double h = 1.0 / n;
    const double a = 0.5 * h * h;

    // Initialize vector b
    Vector b(n + 1);
    for (int i = 0; i <= n; ++i) {
        b(i) = std::cos(M_PI * i * h);
    }

    // Initialize vector u
    Vector u(n + 1); // Assuming the Vector class has a constructor that takes size
    for (int i = 0; i <= n; ++i) {
        u(i) = 0.0; // Initialize with zeros
    }

    // Initialize A_dense
    DenseRowMatrix A_dense(n + 1, n + 1);
    double scale = 1.0 / (h * h);
    for (int i = 0; i < n + 1; ++i) {
        if (i > 0) {
            A_dense(i, i - 1) = -scale;
        }
        A_dense(i, i) = 2.0 * scale;
        if (i < n) {
            A_dense(i, i + 1) = -scale;
        }
    }

    A_dense(0,0)=1/h/h;
    A_dense(n,n)=1/h/h;

    // Initialize A_sparse using the provided code
    int nnz = 3 * (n - 1) + 4;
    int *row_indices = new int[n + 2];
    int *col_indices = new int[nnz];
    double *values = new double[nnz];

    //initialize values
    for (int i = 1; i < n; ++i) {
        values[3 * i] = 2 / h / h;
        values[3 * i - 1] = -1 / h / h;
        values[3 * i + 1] = -1 / h / h;
    }
    values[0] = values[3 * n] = 1 / h / h;
    values[1] = values[3 * n - 1] = -1 / h / h;

    //initialize col_indicies
    for (int i = 0; i < n; i++) {
        col_indices[3 * i + 1] = i + 1;
    }

    for (int i = 0; i <= n; i++) {
        col_indices[3 * i] = i;
    }
    for (int i = 1; i <= n; i++) {
        col_indices[3 * i - 1] = i - 1;
    }
    col_indices[0] = 0;

    //initialize row_indices
    for (int i = 1; i <= n; i++) {
        row_indices[i] = 3 * i - 1;
    }
    row_indices[0] = 0;
    row_indices[n + 1] = 3 * n + 1;

    // Create the SparseMatrixCSR object
    SparseMatrixCSR A_sparse(n + 1, n + 1, nnz, values, col_indices, row_indices);

    // Iteration for Dense Matrix
    int iteration_dense = 0;
    double norm_r_dense = 0.0;
    auto start_dense = std::chrono::high_resolution_clock::now();
    Vector r_dense = b - A_dense * u;
    do {
        r_dense = b - A_dense * u;
        norm_r_dense = norm(r_dense);
        u = u + a * r_dense;
        iteration_dense++;
    } while (norm_r_dense > 1e-3);

    auto end_dense = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_dense = end_dense - start_dense;
    std::cout << "Dense Matrix: Iterations = " << iteration_dense << std::endl;
    std::cout << "Dense Matrix: Time per iteration = " 
              << (elapsed_dense.count() / iteration_dense) << " seconds" << std::endl;
    std::cout << "Final value of norm(r) for Dense Matrix: " << norm_r_dense << std::endl;

    // Reset vector u for Sparse Matrix
    for (int i = 0; i <= n; ++i) {
        u(i) = 0.0; // Reinitialize with zeros
    }

    // Iteration for Sparse Matrix
    int iteration_sparse = 0;
    double norm_r_sparse = 0.0;
    auto start_sparse = std::chrono::high_resolution_clock::now();
    Vector r_sparse = b -  operator*(static_cast<const SparseMatrixCSR &>(A_sparse), u);
    do {
         r_sparse = b -  operator*(static_cast<const SparseMatrixCSR &>(A_sparse), u);
        norm_r_sparse = norm(r_sparse);
        u = u + a * r_sparse;
        iteration_sparse++;
    } while (norm_r_sparse > 1e-3 );

    auto end_sparse = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_sparse = end_sparse - start_sparse;
    std::cout << "Sparse Matrix: Iterations = " << iteration_sparse << std::endl;
    std::cout << "Sparse Matrix: Time per iteration = " 
              << (elapsed_sparse.count() / iteration_sparse) << " seconds" << std::endl;
    std::cout << "Final value of norm(r) for Sparse Matrix: " << norm_r_sparse << std::endl;



        return 0;
}


