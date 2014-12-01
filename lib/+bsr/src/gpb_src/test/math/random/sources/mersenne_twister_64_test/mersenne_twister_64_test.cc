/*
 * Test 64-bit Mersenne Twister random number generator.
 */
#include "lang/array.hh"
#include "lang/exceptions/exception.hh"
#include "lang/exceptions/ex_index_out_of_bounds.hh"
#include "math/random/sources/mersenne_twister_64.hh"

#include <iostream>

using lang::array;
using lang::exceptions::exception;
using lang::exceptions::ex_index_out_of_bounds;
using math::random::sources::mersenne_twister_64;

int main() {
   try {
      array<unsigned long long> seed(0);
      mersenne_twister_64 r(seed);
      for (unsigned long n = 0; n < 1000; n++) {
         std::cout << r.gen_uniform_unsigned_long_long() << " ";
         if (((n+1) % 5) == 0)
            std::cout << "\n";
      }
      for (unsigned long n = 0; n < 1000; n++) {
         std::cout << r.gen_uniform_half_open_double() << " ";
         if (((n+1) % 5) == 0)
            std::cout << "\n";
      }
   } catch (ex_index_out_of_bounds& e) {
      std::cout << e << "\n";
      std::cout << e.index() << "\n";
   } catch (exception& e) {
      std::cout << e << "\n";
   }
   return 0;
}
