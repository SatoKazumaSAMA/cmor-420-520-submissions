#include "abstract_matrix.hpp"

class DenseRowMatrix : public AbstractMatrix {
public:
    DenseRowMatrix(int m, int n);
    ~DenseRowMatrix();

    double operator()(int row, int column) const override;
    double &operator()(int row, int column) override;

private:
    double *_data;
};
