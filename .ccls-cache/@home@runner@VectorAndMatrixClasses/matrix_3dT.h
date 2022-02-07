//============================================================
// FILE: matrix_3dT.h
//============================================================
#ifndef __matrix3d_T_H__
#define __matrix3d_T_H__
#include <cstring>
#include "vector_3dT.h"
template <typename T> class matrix3d;
template <typename T> std::ostream& operator<<(std::ostream& os, const matrix3d<T>& m);\
typedef matrix3d<double> matrix3dD;
typedef matrix3d<float> matrix3dF;
typedef matrix3d<int> matrix3dI;
typedef matrix3d<long> matrix3dL;
template <typename T>
class matrix3d {
public:
 matrix3d();
 matrix3d(const std::string& name, int dims);
 matrix3d(const std::string& name, int dims, const std::initializer_list<vector3d<T>>& li);
 matrix3d(const std::string& name, int dims, const std::initializer_list<T>& li);
//=======================================================================
 matrix3d<T>& operator=(T array[9]);
 matrix3d<T>& operator=(T k);
//=======================================================================
// indexing ops...
 vector3d<T> operator[](int i) const;
 vector3d<T>& operator[](int i);
 T operator()(int row, int col) const; 
 T& operator()(int row, int col);
 T* opengl_memory();
//=======================================================================
 void name(const std::string& name);
 const std::string& name() const;
//============================ LINEAR ALGEBRA =========================
 matrix3d<T>& operator+=(T k);
 matrix3d<T>& operator-=(T k);
 matrix3d<T>& operator*=(T k);
 matrix3d<T>& operator/=(T k);
//=======================================================================
 matrix3d<T>& operator+=(const matrix3d<T>& b); 
 matrix3d<T>& operator-=(const matrix3d<T>& b);
//=======================================================================
 matrix3d<T> operator-();
 matrix3d<T> operator+(const matrix3d<T>& b);
 matrix3d<T> operator-(const matrix3d<T>& b);
//=======================================================================
 friend matrix3d operator+(const matrix3d& a, T k) {
return matrix3d(std::to_string(k) + "+" + a.name(), 3,
{ a[0] + k, a[1] + k, a[2] + k});
 }
 friend matrix3d operator+(T k, const matrix3d& a) { return a + k; }
 friend matrix3d operator-(const matrix3d& a, T k) { return a + -k; }
 friend matrix3d operator-(T k, const matrix3d& a) {
// implement code here
 }
 friend matrix3d operator*(const matrix3d& a, T k) {
// implement code here
 }
 friend matrix3d<T> operator*(T k, const matrix3d& a) { return a * k; }
 friend matrix3d operator/(const matrix3d& a, T k) {
// implement code here
 }
//=======================================================================
 friend matrix3d operator*(const matrix3d& m, const vector3d<T>& v) {
// implement code here
 }
 friend matrix3d operator*(const vector3d<T>& v, const matrix3d& m) {
// implement code here
 }
 matrix3d<T> operator*(const matrix3d<T>& b);
//=======================================================================
 matrix3d<T> transpose() const; // create a new matrix transpose()
 T determinant() const;
 T trace() const;
//=======================================================================
 matrix3d<T> minors() const; // see defn
 matrix3d<T> cofactor() const; // (-1)^(i+j)*minors()(i, j)
 matrix3d<T> adjugate() const; // cofactor.transpose()
 matrix3d<T> inverse() const; // adjugate()/determinant()
//=======================================================================
 static matrix3d<T> identity(int dims); // identity matrix
 static matrix3d<T> zero(int dims); // zero matrix
//=======================================================================
 bool operator==(const matrix3d<T>& b) const;
 bool operator!=(const matrix3d<T>& b) const;
//=======================================================================
 friend std::ostream& operator<< <> (std::ostream& os, const matrix3d<T>& m);
private: 
 void check_equal_dims(const matrix3d<T>& v) const;
 void check_bounds(int i) const;
 void swap(T& x, T& y);
private: 
 std::string name_;
 int dims_;
 vector3d<T> cols_[4];
 T data_[16];
};
//=================================================================================================
template <typename T> matrix3d<T>::matrix3d() : matrix3d("", 3) {} // 3d default dims
template <typename T> matrix3d<T>::matrix3d(const std::string& name, int dims)
: name_(name), dims_(dims) {
 for (int i = 0; i < 4; ++i) { cols_[i].name("col" + std::to_string(i)); }
 std::memset(data_, 0, 16 * sizeof(T));
}
template <typename T> matrix3d<T>::matrix3d(const std::string& name, int dims,
const std::initializer_list<vector3d<T>>& li)
: matrix3d(name, dims) {
 int i = 0;
 for (vector3d<T> value : li) {
if (i > dims_) { break; }
cols_[i++] = value;
 }
}
template <typename T> matrix3d<T>::matrix3d(const std::string& name, int dims,
const std::initializer_list<T>& li)
: matrix3d(name, dims) {
 int i = 0;
 for (T value : li) {
cols_[i/3][i % 3] = value;
	 ++i;
 }
}
//=================================================================================================
template <typename T> matrix3d<T>& matrix3d<T>::operator=(T array[9]) {
 for (int i = 0; i < 3; ++i) {
for (int j = 0; j < 3; ++i) {
cols_[i][j] = array[i + j];
}
 }
 return *this;
}
template <typename T> matrix3d<T>& matrix3d<T>::operator=(T k) {
 for (int i = 0; i < 3; ++i) {
for (int j = 0; j < 3; ++j) {
cols_[i][j] = k;
}
 }
 return *this;
}
//=================================================================================================
template <typename T> vector3d<T> matrix3d<T>::operator[](int i) const {
 check_bounds(i); return cols_[i];
}
template <typename T> vector3d<T>& matrix3d<T>::operator[](int i) {
 check_bounds(i); return cols_[i];
}
template <typename T> T matrix3d<T>::operator()(int row, int col) const {
// implement code here
}
template <typename T> T& matrix3d<T>::operator()(int row, int col) {
// implement code here
}
template <typename T> T* matrix3d<T>::opengl_memory() { // constant ptr
// implement code here
}
//=================================================================================================
template <typename T> void matrix3d<T>::name(const std::string& name) { name_ = name; }
template <typename T> const std::string& matrix3d<T>::name() const { return name_; }
//=================================== LINEAR ALGEBRA ================================
template <typename T> matrix3d<T>& matrix3d<T>::operator+=(T k) {
 const matrix3d<T>& a = *this;
 name_ = std::to_string(k) + "+" + name_;
 for (int i = 0; i < 4; ++i) { a[i] += k; }
 return *this;
}
template <typename T> matrix3d<T>& matrix3d<T>::operator-=(T k) { *this += -k; return *this; }
template <typename T> matrix3d<T>& matrix3d<T>::operator*=(T k) {
// implement code here
}
template <typename T> matrix3d<T>& matrix3d<T>::operator/=(T k) {
// implement code here
}
//=================================================================================================
template <typename T> matrix3d<T>& matrix3d<T>::operator+=(const matrix3d<T>& b) {
// implement code here
}
template <typename T> matrix3d<T>& matrix3d<T>::operator-=(const matrix3d<T>& b) {
// implement code here
}
//=================================================================================================
template <typename T> matrix3d<T> matrix3d<T>::operator-() {
 const matrix3d<T>& a = *this;
	 return matrix3d<T>("-" + name_, 3, {-a[0], -a[1], -a[2]});
}
template <typename T> matrix3d<T> matrix3d<T>::operator+(const matrix3d<T>& b) {
 const matrix3d<T>& a = *this;
 check_equal_dims(b);
 return matrix3d<T>(name_ + "+" + b.name_, dims_, {a[0] + b[0], a[1] + b[1], a[2] + b[2]});
}
template <typename T> matrix3d<T> matrix3d<T>::operator-(const matrix3d<T>& b) {
// implement code here
}
//=================================================================================================
template <typename T> matrix3d<T> matrix3d<T>::operator*(const matrix3d<T>& b) {
 const matrix3d<T>& a = *this;
 return matrix3d<T>(a.name_ + "*" + b.name_, 3, {
a(0,0)*b(0,0) + a(0,1)*b(1,0) + a(0,2)*b(2,0),
a(1,0)*b(0,0) + a(1,1)*b(1,0) + a(1,2)*b(2,0),
a(2,0)*b(0,0) + a(2,1)*b(1,0) + a(2,2)*b(2,0),
a(0,0)*b(0,1) + a(0,1)*b(1,1) + a(0,2)*b(2,1),
a(1,0)*b(0,1) + a(1,1)*b(1,1) + a(1,2)*b(2,1),
a(2,0)*b(0,1) + a(2,1)*b(1,1) + a(2,2)*b(2,1),
a(0,0)*b(0,2) + a(0,1)*b(1,2) + a(0,2)*b(2,2),
a(1,0)*b(0,2) + a(1,1)*b(1,2) + a(1,2)*b(2,2),
a(2,0)*b(0,2) + a(2,1)*b(1,2) + a(2,2)*b(2,2)} );
}
//=================================================================================================
template <typename T> matrix3d<T> matrix3d<T>::transpose() const {
 const matrix3d<T>& m = *this;
// implement code here
}
template <typename T> T matrix3d<T>::determinant() const {
// implement code here
}
template <typename T> T matrix3d<T>::trace() const {
 const matrix3d<T>& m = *this;
 return m(0,0) + m(1,1) + m(2,2);
}
//=================================================================================================
// | | e f | | d f | | d e | | Matrix of minors
// | | h i | | g i | | g h | |
// | |
// | | b c | | a c | | a b | |
// | | h i | | g i | | g h | |
// | |
// | | b c | | a c | | a b | |
// | | e f | | d f | | d e | |
// ||
//----------------------------------------------------------------
template <typename T> matrix3d<T> matrix3d<T>::minors() const {
 const matrix3d<T>& m = *this;
 return matrix3d<T>("Min(" + name_ + ")", 3, {
(m(1,1)*m(2,2) - m(1,2)*m(2,1)),
(m(0,1)*m(2,2) - m(0,2)*m(2,1)),
(m(0,1)*m(1,2) - m(0,2)*m(1,1)),
(m(1,0)*m(2,2) - m(1,2)*m(2,0)),
(m(0,0)*m(2,2) - m(0,2)*m(2,0)),
(m(0,0)*m(1,2) - m(0,2)*m(1,0)),
(m(1,0)*m(2,1) - m(1,1)*m(2,0)),
(m(0,0)*m(2,1) - m(0,1)*m(2,0)),
(m(0,0)*m(1,1) - m(0,1)*m(1,0)) });
}
template <typename T> matrix3d<T> matrix3d<T>::cofactor() const {
// implement code here
}
template <typename T> matrix3d<T> matrix3d<T>::adjugate() const {
// implement code here
}
template <typename T> matrix3d<T> matrix3d<T>::inverse() const {
// implement code here
}
//=================================================================================================
template <typename T> matrix3d<T> matrix3d<T>::identity(int dims) {
// implement code here
}
template <typename T> matrix3d<T> matrix3d<T>::zero(int dims) {
// implement code here
}
template <typename T> bool matrix3d<T>::operator==(const matrix3d<T>& b) const {
 check_equal_dims(b);
 const matrix3d<T>& a = *this;
 return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}
template <typename T> bool matrix3d<T>::operator!=(const matrix3d<T>& b) const {
 return !(*this == b);
}
//=================================================================================================
template <typename T> std::ostream& operator<<(std::ostream& os, const matrix3d<T>& m) {
 os << "<'" << m.name_ << "', ";
 for (int i = 0; i < 3; ++i) { os << m.cols_[i]; }
 os << "> OR by rows...\n";
 for (int i = 0; i < 3; ++i) {
for (int j = 0; j < 3; ++j) {
os << m(i, j) << " ";
}
os << "\n";
 }
 return os << ">";
}
//=================================================================================================
template <typename T> void matrix3d<T>::check_equal_dims(const matrix3d<T>& v) const {
 if (dims_ != v.dims_) { throw new std::invalid_argument("matrix3d dims mismatch"); }
}
template <typename T> void matrix3d<T>::check_bounds(int i) const {
 if (i > dims_) {
throw new std::invalid_argument("out of bounds");
 }
}
template <typename T> void matrix3d<T>::swap(T& x, T& y) {
 T temp = x; x = y; y = temp;
}
//============================================================
// end of file: matrix_3dT.h
//============================================================
#endif