#include "abstract_matrix.hpp"

Vector operator*(AbstractMatrix & A, Vector & x){
  Vector out(x.length());
  for (int i = 0; i < A.num_rows(); ++i){
    out(i) = 0.0;
    for (int j = 0; j < A.num_rows(); ++j){
      out(i) += A(i,j) * x(j);
    }
  }
  return out;
}
