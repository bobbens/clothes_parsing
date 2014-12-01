/*
 * Queue set test.
 */
#include "collections/queue_set.hh"
#include "functors/comparable_functors.hh"
#include "lang/array.hh"
#include "lang/exceptions/exception.hh"
#include "math/random/generators/rand_gen_uniform.hh"

#include <iostream>

using collections::queue_set;
using functors::compare_functors;
using lang::array;
using lang::exceptions::exception;
using math::random::generators::rand_gen_uniform;

int main() {
   try {
      unsigned long N = 100;
      queue_set<double> q(
         compare_functors<double>::f_compare(),
         compare_functors<double>::f_compare_memloc()
      );
      rand_gen_uniform<> r;
      array<double> a(N);
      for (unsigned long n = 0; n < N; n++)
         a[n] = r.generate();
      for (unsigned long n = 0; n < N; n++) {
         unsigned long ind = static_cast<unsigned long>(r.generate() * (N-1));
         q.add(a[ind]);
      }
      for (unsigned long n = 0; n < (N/2); n++) {
         unsigned long ind = static_cast<unsigned long>(r.generate() * (N-1));
         q.remove(a[ind]);
      }
      if (!q.is_empty()) {
         double x = q.dequeue();
         while (!q.is_empty()) {
            double y = q.dequeue();
            if (x > y) {
               std::cout << "failed\n";
               return 0;
            }
            x = y;
         }
      }
      std::cout << "passed\n";
   } catch (exception& e) {
      std::cout << e << "\n";
      throw;
   }
   return 0;
}
