/*
 * Signal processing libarary test.
 */
#include "lang/exceptions/exception.hh"
#include "math/libraries/lib_signal.hh"
#include "math/complex.hh"
#include "math/math.hh"
#include "math/matrices/cmatrix.hh"
#include "math/matrices/matrix.hh"
#include "io/streams/cout.hh"

using lang::exceptions::exception;
using math::complex;
using math::libraries::lib_signal;
using math::matrices::cmatrix;
using math::matrices::matrix;
using io::streams::cout;

int main() {
   try {
      matrix<>   m = matrix<>::ramp(1,8);
      cmatrix<> cm = m + complex<>(0,1)*matrix<>::ramp(-0.3,0.1,0.41);
      m.transpose();
      cm.transpose();
      cout << m << "\n";
      cout << cm << "\n";
      cout << "--- dft/idft ---\n";
      cout << lib_signal::dft(m,1) << "\n";
      cout << lib_signal::idft(lib_signal::dft(m,1),1) << "\n";
      cout << "--- fft/ifft ---\n";
      cout << lib_signal::fft(m,1) << "\n";
      cout << lib_signal::ifft(lib_signal::fft(m,1),1) << "\n";
      cout << "--- dct/idct (real) ---\n";
      cout << lib_signal::dct(m,1) << "\n";
      cout << lib_signal::idct(lib_signal::dct(m,1),1) << "\n";
      cout << "--- dct/idct (complex) ---\n";
      cout << lib_signal::dct(cm,1) << "\n";
      cout << lib_signal::idct(lib_signal::dct(cm,1),1) << "\n";
      cout << "--- dst/idst (real) ---\n";
      cout << lib_signal::dst(m,1) << "\n";
      cout << lib_signal::idst(lib_signal::dst(m,1),1) << "\n";
      cout << "--- dst/idst (complex) ---\n";
      cout << lib_signal::dst(cm,1) << "\n";
      cout << lib_signal::idst(lib_signal::dst(cm,1),1) << "\n";
      cout << "--- hilbert ---\n";
      cout << lib_signal::hilbert(m,1) << "\n";
   } catch (exception& e) {
      cout << e << "\n";
   }
   return 0;
}
