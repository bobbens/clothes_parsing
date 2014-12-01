/*
 * Test complex numbers.
 * Right now, this really just tests parsing of complex.hh
 */
#include "math/complex.hh"
#include <iostream>

using math::complex;

int main() {
   complex<> z;
   complex<> z0(1, 2);
   complex<> z1(5, 1);
   std::cout << z0 + z1 << "\n";
   std::cout << z.real() << "\n";
   return 0;
}
