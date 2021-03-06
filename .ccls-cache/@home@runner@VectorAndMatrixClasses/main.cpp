//============================================================
// CPSC484 - Spring 2022, CSUF, Dr. William McCarthy
// Student: Pham, Trong
// Student: Nguyen, Michael
// file: main.cpp
//============================================================
#include <iostream>
#include <cstring>
#include <initializer_list>
#include <cassert>
#include "matrix_3dT.h"
#include "vector_3dT.h"
#include "quaternion_T.h"
template <typename T> void print(T v);
template <typename T> void show_vect(T v);
template <typename T> void show_mat(T m);
void test_vectors();
void test_matrices();
void test_matrices_and_vectors();
void test_quaternions();
int main(int argc, const char * argv[]) {
 test_vectors();
 test_matrices();
 test_matrices_and_vectors();
 print("... program completed ! ...\n");
 return 0;
}
void test_quaternions() {
    print("\n====================  TESTING QUATERNIONS  ========================");
    quaternion<double>::run_tests();
    print("...test_matrices_and_vectors assertions passed");
    print("====================  FINISHED testing quaternions  ========================");
}

void test_matrices_and_vectors() {
 print("\n==================== TESTING MATRICES and VECTORS ========================");
 //vector3dD flushing("fp", 2, {1,5});//for flushing purpose only
 vector3dD p("p",2,{1,2});
 show_vect(p);
	
 matrix3dD m("m", 2, {1, 2, 3, 4});
 show_mat(m);
 
 assert(p * m == m * p);
	
 vector3dD q("q", 3, {1, 2, 3});
 matrix3dD n("n", 3, {1, 2, 3, 4, 5, 6, 7, 8, 9});
 show_vect(q);
 show_mat(n);
 assert(q * n == n * q);
 print("...test_matrices_and_vectors assertions passed");
 print("==================== FINISHED testing matrices and vectors ========================");
}
void test_vectors() {
 print("\n==================== TESTING VECTORS ========================");
 
 vector3dD u("u", 3, {1, 2, 4});
 vector3dD v("v", 3, {8, 16, 32});
 show_vect(u);
 show_vect(v);
 
 vector3dD i("i", 3, {1, 0, 0}), j("j", 3, {0, 1, 0}), k("k", 3, {0, 0, 1});
 show_vect(i);
 show_vect(j);
 show_vect(k);
 
 vector3dD w(3 * i + 4 * j - 2 * k);
 show_vect(w);
 
 assert(u == u);
 assert(u != v);

 
 assert(u + v == v + u);
 assert(u - v == -(v - u));
 
 assert(-(-u) == u);
 assert(3.0 + u == u + 3.0);
 assert(3.0 * u == u * 3.0);
 assert((u - 3.0) == -(3.0 - u));
 
 assert((5.0 * u) / 5.0 == u);
 
 assert(u + vector3dD::zero() == u);
 
 assert(i.dot(j) == j.dot(k) == k.dot(i) == 0);
 assert(i.cross(j) == k);
 assert(j.cross(k) == i);
 assert(k.cross(i) == j);
 assert(u.cross(v) == -v.cross(u));
 
 assert(u.cross(v + w) == u.cross(v) + u.cross(w));
 assert((u.cross(v)).dot(u) == 0);
 print(i.angle(j));
 print(M_PI/2);
 
 assert(i.angle(j) == M_PI_2);
 assert(j.angle(k) == M_PI_2);
 assert(k.angle(i) == M_PI_2);
 vector3dD uhat = u / u.magnitude(); // unit vector in u direction
 show_vect(u);
 show_vect(uhat);
 print(uhat.magnitude());
 
 assert(uhat.magnitude() - 1.0 < 1.0e-10);
 print("...test vectors assertions passed");
 print("==================== FINISHED testing vectors ========================");
 
}
void test_matrices() {
 print("\n==================== TESTING MATRICES ========================");
 matrix3dD a("a", 3, {3, 2, 0, 0, 0, 1, 2, -2, 1});
 print(a);
		
 matrix3dD b("b", 3, {1, 0, 5, 2, 1, 6, 3, 4, 0});
 print(b);
 
 matrix3dD ainv = a.inverse();
 print(ainv);
 
 matrix3dD binv = b.inverse();
 print(binv);
 
 std::cout<<"a*ainv"<<std::endl;
 print(a * ainv);
 print(b * binv);
 
 assert( a * ainv == matrix3dD::identity(3));
	
 assert(a * ainv == ainv * a);
 assert(b * binv == matrix3dD::identity(3));
 assert(b * binv == binv * b);

 assert(a.transpose().transpose() == a);
 assert(a.transpose().determinant() == a.determinant());
 assert(a + b == b + a);
 assert(a - b == -(b - a));
 assert(3.0 + a == a + 3.0);
 assert(3.0 * a == a * 3.0);
 assert((a + 3.0) - 3.0 == a);
 assert((3.0 * a) / 3.0 == a);
	
 assert(-(-a) == a);
 matrix3dD zerod("zerod", 3, {1, 2, 3, 4, 5, 6, 7, 8, 9});
 assert(zerod.determinant() == 0);
 print("...test matrices assertions passed");
 print("==================== FINISHED testing matrices ========================");
}
template <typename T> void print(T v) {
 std::cout << v << std::endl;
}
template <typename T> void show_vect(T v) {
 std::cout << v.name() << " is: " << v << std::endl;
}
template <typename T> void show_mat(T m) {
 std::cout << m.name() << " is: " << m << std::endl;
}
