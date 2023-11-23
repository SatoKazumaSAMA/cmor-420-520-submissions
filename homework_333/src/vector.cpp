#include "vector.hpp"
#include <cmath>
Vector & Vector::operator+=(Vector x){
  for (int i = 0; i < x.length(); ++i){
    _data[i] += x(i);
  }
  return *this; // pointer to the instance 
}

Vector & Vector::operator=(const Vector & copy_from){
  if (_length != copy_from.length()){
    std::cout << "ERROR: vectors are not the same length";
    std::cout << std::endl;
  }  
  for (int i = 0; i < _length; ++i){
    _data[i] = copy_from(i);
  }

  return *this;
}

Vector::Vector(const Vector & copy_from){
  _length = copy_from.length();
  _data = new double [_length];
  for (int i = 0; i < _length; ++i){
    _data[i] = copy_from(i);
  }
}

void Vector::print(std::string variable_name){
  std::cout << variable_name << " = " << std::endl;
  for (int i = 0; i < _length; ++i){
    std::cout << _data[i] << std::endl;
  }
}

Vector operator+(const Vector & x, const Vector & y){
  Vector out(x.length());
  for (int i = 0; i < out.length(); ++i){
    out(i) = x(i) + y(i);
  }
  return out;
}

Vector operator-(const Vector & x, const Vector & y){
  return (x + -1.0 * y);
}

Vector operator*(double x, const Vector & y){
  Vector out(y.length());
  for (int i = 0; i < out.length(); ++i){
    out(i) = x * y(i);
  }
  return out;  
}





Vector operator*(const Vector & x, double y){
  return y * x;
}



