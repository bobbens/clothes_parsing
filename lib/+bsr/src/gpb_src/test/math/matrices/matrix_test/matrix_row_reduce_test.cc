/*
 * Matrix test.
 */
#include "io/streams/cout.hh"
#include "lang/array.hh"
#include "math/matrices/exceptions/ex_matrix_singular.hh"
#include "math/matrices/matrix.hh"

#include "lang/exceptions/ex_index_out_of_bounds.hh"

using io::streams::cout;
using lang::array;
using math::matrices::exceptions::ex_matrix_singular;
using math::matrices::matrix;

int main() {
   try {
      matrix<> m(4,4);
      m(0,0) = 0.7468; m(0,1) = 0.4186; m(0,2) = 0.6721; m(0,3) = 0.3795;
      m(1,0) = 0.4451; m(1,1) = 0.8462; m(1,2) = 0.8381; m(1,3) = 0.8318;
      m(2,0) = 0.9318; m(2,1) = 0.5252; m(2,2) = 0.0196; m(2,3) = 0.5028;
      m(3,0) = 0.4660; m(3,1) = 0.2026; m(3,2) = 0.6813; m(3,3) = 0.7095;
      matrix<> a(m);
      for (unsigned long n = 0; n < 4; n++)
         a(2,n) = a(1,n);
      cout << "--- input matrices ---\n";
      cout << "m = " << m << "\n";
      cout << "a = " << a << "\n";
      cout << "--- ref, rref ---\n";
      cout << "ref(m) = " << ref(m) << "\n";
      cout << "rref(m) = " << rref(m) << "\n";
      cout << "ref(a) = " << ref(a) << "\n";
      cout << "rref(a) = " << rref(a) << "\n";
      cout << "--- det, rank ---\n";
      cout << "m.det()  = " << m.det() << "\n";
      cout << "m.rank() = " << m.rank() << "\n";
      cout << "a.det()  = " << a.det() << "\n";
      cout << "a.rank() = " << a.rank() << "\n";
      cout << "--- inv ---\n";
      matrix<> m_inv = inv(m);
      cout << "inv(m) = " << m_inv << "\n";
      cout << "inv(inv(m)) = " << inv(m_inv) << "\n";
      cout << "m*inv(m) = " << m * m_inv << "\n";
      cout << "inv(m)*m = " << m_inv * m << "\n";
      cout << "inv(m).det() = " << m_inv.det() << "\n";
      cout << "inv(m).rank() = " << m_inv.rank() << "\n";
      try {
         cout << "inv(a) = " << inv(a) << "\n";
      } catch (math::matrices::exceptions::ex_matrix_singular& e) {
         cout << "a is singular\n";
      }
   } catch (lang::exceptions::ex_index_out_of_bounds& e) {
      cout << e << "\n";
      cout << e.index() << "\n";
   }
   return 0;
}
