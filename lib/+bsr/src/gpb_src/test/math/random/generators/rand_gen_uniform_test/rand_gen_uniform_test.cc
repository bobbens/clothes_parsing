/*
 * Test uniform random generator.
 */
#include "lang/exceptions/exception.hh"
#include "lang/exceptions/ex_index_out_of_bounds.hh"
#include "math/random/generators/rand_gen_uniform.hh"

#include <iostream>

using lang::exceptions::exception;
using lang::exceptions::ex_index_out_of_bounds;
using math::random::generators::rand_gen_uniform;

int main(int argc, char** argv) {
   try {
      for (unsigned long m = 0; m < 10; m++) {
         rand_gen_uniform<> r;
         for (unsigned long n = 0; n < 1000; n++) {
            if (argc > 1) {
               std::cout << r.generate() << " ";
               if (((n+1) % 5) == 0)
                  std::cout << "\n";
            }
         }
      }
   } catch (ex_index_out_of_bounds& e) {
      std::cout << e << "\n";
      std::cout << e.index() << "\n";
   } catch (exception& e) {
      std::cout << e << "\n";
   }
   std::cout << "rand_gen_uniform_test completed\n";
   return 0;
}
