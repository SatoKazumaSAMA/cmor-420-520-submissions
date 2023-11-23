#include<iostream>

#ifndef _VECTOR
#define _VECTOR

class Vector{

public:
  // constructor
  Vector(int length){
    _length = length;
    _data = new double [length];
  }

  // destructor  
  ~Vector(){ delete [] _data; }

  // copy constructor
  Vector(const Vector & copy_from);

  // copy move constructor
  Vector & operator=(const Vector & copy_from);

  int length() const{ return _length; };

  void print(std::string variable_name);

  double & operator()(const int index) const {  
    return _data[index];
  }
  
  Vector & operator+=(Vector x);

  double norm() const;

  Vector & mul_nalloc(double x);
  Vector & add_nalloc(Vector &x);
  Vector& minus_equal_nalloc(const Vector& x,const Vector &b);
  
private:  

  int _length;
  double * _data;
  
};

#endif

Vector operator+(const Vector & x, const Vector & y);
Vector operator-(const Vector & x, const Vector & y);
Vector operator*(double x, const Vector & y);
Vector operator*(const Vector & x, double y);
