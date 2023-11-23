#include "vector.hpp"

#ifndef _ABSTRACTMATRIX
#define _ABSTRACTMATRIX

class AbstractMatrix{
public:
  AbstractMatrix(int m, int n){
    _rows = m;
    _columns = n;
  }
  int num_rows()const{ return _rows; };
  int num_columns()const{ return _columns; };

  virtual double operator()(int row, int column) const{
    throw std::logic_error("`double operator()` not implemented"
			   "for AbstractMatrix");
  }
  
  virtual double & operator()(int row, int column){
    throw std::logic_error("`double & operator()` not implemented"
			   "for AbstractMatrix");
  }

private:
  int _rows;
  int _columns;
};

#endif

Vector operator*(AbstractMatrix & A, Vector & x);



