/*
 * Matrix test.
 */
#include "io/streams/cout.hh"
#include "lang/array.hh"
#include "math/matrices/matrix.hh"

#include "lang/exceptions/ex_index_out_of_bounds.hh"

using io::streams::cout;
using lang::array;
using math::matrices::matrix;

int main() {
try {
   matrix<> m(3,2);
   m(0,0) = 0.1;
   m(0,1) = 0.02;
   m(1,0) = 1.0;
   m(1,1) = 2.0;
   m(2,0) = 10;
   m(2,1) = 30;
   matrix<> m_ramp = matrix<>::ramp(-0.5, 0.25, 0.5);
   array<unsigned long> dims1(1);
   dims1[0] = m_ramp.size();
   matrix<> m_ramp1 = reshape(m_ramp,dims1);
   cout << m << "\n";
   cout << m_ramp << "\n";
   cout << m_ramp1 << "\n";
   cout << "--- convolution 1D (full) ---\n";
   cout << conv(m_ramp1,m_ramp1) << "\n";
   cout << "--- convolution 1D (cropped) ---\n";
   cout << conv_crop(m_ramp1,m_ramp1) << "\n";
   cout << "--- convolution 1D (cropped strict) ---\n";
   cout << conv_crop_strict(m_ramp1,m_ramp1) << "\n";
   cout << "--- convolution 2D (full) ---\n";
   cout << conv(m,m) << "\n";
   cout << conv(m,transpose(m)) << "\n";
   cout << conv(m,m_ramp) << "\n";
   cout << conv(m_ramp,m) << "\n";
   cout << conv(m_ramp,transpose(m_ramp)) << "\n";
   cout << "--- convolution 2D (cropped) ---\n";
   cout << conv_crop(m,m) << "\n";
   cout << conv_crop(m,transpose(m)) << "\n";
   cout << conv_crop(m,m_ramp) << "\n";
   cout << conv_crop(m_ramp,m) << "\n";
   cout << conv_crop(m_ramp,transpose(m_ramp)) << "\n";
   cout << "--- convolution 2D (cropped strict) ---\n";
   cout << conv_crop_strict(m,m) << "\n";
   cout << conv_crop_strict(m,transpose(m)) << "\n";
   cout << conv_crop_strict(m,m_ramp) << "\n";
   cout << conv_crop_strict(m_ramp,m) << "\n";
   cout << conv_crop_strict(m_ramp,transpose(m_ramp)) << "\n";
   cout << "--- min, max, sum, prod ---\n";
   cout << min(m,0) << "\n";
   cout << min(m,1) << "\n";
   cout << max(m,0) << "\n";
   cout << max(m,1) << "\n";
   cout << sum(m,0) << "\n";
   cout << sum(m,1) << "\n";
   cout << prod(m,0) << "\n";
   cout << prod(m,1) << "\n";
   cout << "--- cumsum, cumprod ---\n";
   cout << cumsum(m,0) << "\n";
   cout << cumsum(m,1) << "\n";
   cout << cumprod(m,0) << "\n";
   cout << cumprod(m,1) << "\n";
   cout << "--- mean, var ---\n";
   cout << mean(m) << "\n";
   cout << var(m) << "\n";
   cout << mean(m,0) << "\n";
   cout << mean(m,1) << "\n";
   cout << var(m,0) << "\n";
   cout << var(m,1) << "\n";
   cout << "--- gradient ---\n";
   cout << gradient(m,0) << "\n";
   cout << gradient(m,1) << "\n";
   cout << "--- reverse ---\n";
   cout << m.reverse(0) << "\n";
   cout << m.reverse(1) << "\n";
   cout << "--- vertcat, horzcat, concat ---\n";
   cout << vertcat(m,m) << "\n";
   cout << horzcat(m,m) << "\n";
   cout << concat(m,m,4) << "\n";
   cout << "--- resize, transpose ---\n";
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
   cout << "--- repmat ---\n";
   cout << repmat(m,order) << "\n";
   cout << "--- sort ---\n";
   cout << m << "\n";
   cout << m.sort_idx(0) << "\n";
   cout << m << "\n";
   cout << m.sort_idx(1) << "\n";
   cout << m << "\n";
} catch (lang::exceptions::ex_index_out_of_bounds& e) {
   cout << e << "\n";
   cout << e.index() << "\n";
}
   return 0;
}
