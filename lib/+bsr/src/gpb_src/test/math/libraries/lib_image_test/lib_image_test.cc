/*
 * Image processing library test.
 */
#include "lang/exceptions/exception.hh"
#include "lang/exceptions/ex_index_out_of_bounds.hh"
#include "math/libraries/lib_image.hh"
#include "math/math.hh"
#include "math/matrices/matrix.hh"
#include "io/streams/cout.hh"

using lang::exceptions::exception;
using lang::exceptions::ex_index_out_of_bounds;
using math::libraries::lib_image;
using math::matrices::matrix;
using io::streams::cout;

int main() {
   try {
      cout << "--- 2D matrix resampling ---\n";
      matrix<> g = lib_image::gaussian_2D();
      cout << g << "\n";
      matrix<> g_big = lib_image::resample_2D(g,2);
      cout << g_big << "\n";
      matrix<> g_resampled = lib_image::resample_2D(g_big,0.5);
      cout << g_resampled << "\n";
      unsigned long n = 2;
      matrix<> g_two = lib_image::resample_2D(g_resampled,n,n);
      cout << g_two << "\n";
      n = 1;
      matrix<> g_one = lib_image::resample_2D(g_resampled,n,n);
      cout << g_one << "\n";
      cout << "--- 2D matrix rotation ---\n";
      cout << lib_image::rotate_2D(g,M_PIl/4) << "\n";
      cout << "--- guassian kernels ---\n";
      cout << lib_image::gaussian_2D(1,2,M_PIl/8,0) << "\n";
      cout << lib_image::gaussian_2D(1,2,M_PIl/8,1) << "\n";
      cout << lib_image::gaussian_2D(1,2,M_PIl/8,2) << "\n";
      cout << lib_image::gaussian_2D(1,2,M_PIl/8,2,true) << "\n";
      cout << "--- nonmax ---\n";
      matrix<> m = lib_image::gaussian_2D(1,2,M_PIl/8,2,true);
      m -= min(m);
      cout << m << "\n";
      cout << lib_image::nonmax(m) << "\n";
      cout << lib_image::nonmax(m,0,true) << "\n";
      cout << lib_image::nonmax(m,1,true) << "\n";
   } catch (ex_index_out_of_bounds& e) {
      cout << e.index() << "\n";
      cout << e << "\n";
   } catch (exception& e) {
      cout << e << "\n";
   }
   return 0;
}
