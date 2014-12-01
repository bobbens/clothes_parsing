/*
 * Complex matrix test.
 */
#include "io/streams/cout.hh"
#include "lang/array.hh"
#include "math/complex.hh"
#include "math/matrices/cmatrix.hh"

#include <iomanip>

#include "lang/exceptions/ex_index_out_of_bounds.hh"

using io::streams::cout;
using lang::array;
using math::complex;
using math::matrices::cmatrix;

int main() {
try {
//   cmatrix<> m = cmatrix<>::random(4,5);
   cmatrix<> m(3,2);
   m(0,0) = complex<>(0.1,0.5);
   m(0,1) = complex<>(0.02,0);
   m(1,0) = complex<>(1.0,-1.0);
   m(1,1) = complex<>(2.0,0);
   m(2,0) = complex<>(10,0);
   m(2,1) = complex<>(30,0);
   cout << m << "\n";
   cout << cumsum(m,0) << "\n";
   cout << cumsum(m,1) << "\n";
   cout << cumprod(m,0) << "\n";
   cout << cumprod(m,1) << "\n";
   cout << m.reverse(0) << "\n";
   cout << m.reverse(1) << "\n";
   cout << vertcat(m,m) << "\n";
   cout << horzcat(m,m) << "\n";
   cout << concat(m,m,4) << "\n";
   cout << resize(m,2,2) << "\n";
   cout << resize(m,4,4) << "\n";
   cout << transpose(resize(m,4,4)) << "\n";
   array<unsigned long> order(3);
   order[0] = 2;
   order[1] = 1;
   order[2] = 0;
   cout << permute_dimensions(resize(m,4,4),order) << "\n";
   order[1] = 3;
   order.resize(2);
   cout << repmat(m,order) << "\n";
} catch (lang::exceptions::ex_index_out_of_bounds& e) {
   cout << e << "\n";
   cout << e.index() << "\n";
}
   return 0;
}
