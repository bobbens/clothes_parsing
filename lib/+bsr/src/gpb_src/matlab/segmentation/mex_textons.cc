/*
 * Textons. 
 */
#include "collections/pointers/auto_collection.hh"
#include "collections/array_list.hh"
#include "io/streams/cout.hh"
#include "lang/array.hh"
#include "lang/exceptions/exception.hh"
#include "lang/exceptions/ex_index_out_of_bounds.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/libraries/lib_image.hh"
#include "math/math.hh"
#include "math/matrices/matrix.hh"

#include "concurrent/threads/thread.hh"

#include <time.h>
#include <mex.h>

using collections::pointers::auto_collection;
using collections::array_list;
using io::streams::cout;
using lang::array;
using lang::exceptions::exception;
using lang::exceptions::ex_index_out_of_bounds;
using lang::pointers::auto_ptr;
using math::libraries::lib_image;
using math::matrices::matrix;

using concurrent::threads::thread;

/********************************** 
 * Matlab matrix conversion routines.
 **********************************/

/*
 * Get a string from an mxArray.
 */
char *mexGetString(const mxArray *arr) {
   char *string;
   int buflen;
   buflen = (mxGetM(arr) * mxGetN(arr) * sizeof(mxChar)) + 1;
   string = new char[buflen];
   mxGetString(arr, string, buflen);
   return string;
}

/*
 * Create a single element Matlab double matrix.
 */
mxArray* mxCreateScalarDouble(double value) {
   mxArray* pa = mxCreateDoubleMatrix(1, 1, mxREAL);
   *mxGetPr(pa) = value;
   return pa;
}

/* 
 * Convert an mxArray to a matrix.
 */
matrix<> to_matrix(const mxArray *a) {
   unsigned long mrows = static_cast<unsigned long>(mxGetM(a));
   unsigned long ncols = static_cast<unsigned long>(mxGetN(a));
   double *data = mxGetPr(a);
   matrix<> m(mrows, ncols);
   for (unsigned long r = 0; r < mrows; r++) {
      for (unsigned long c = 0; c < ncols; c++) {
         m(r,c) = data[(c*mrows) + r];
      }
   }
   return m;
}
   
/*
 * Convert a 2D matrix to an mxArray.
 */
mxArray* to_mxArray(const matrix<>& m) {
   unsigned long mrows = m.size(0);
   unsigned long ncols = m.size(1);
   mxArray *a = mxCreateDoubleMatrix(
      static_cast<int>(mrows),
      static_cast<int>(ncols),
      mxREAL
   );
   double *data = mxGetPr(a);
   for (unsigned long r = 0; r < mrows; r++) {
      for (unsigned long c = 0; c < ncols; c++) {
         data[(c*mrows) + r] = m(r,c);
      }
   }
   return a;
}

/*
 * Matlab interface.
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
   try {
      /* set filter parameters */
      unsigned long n_ori = 8;
      double sigma = 4;
      unsigned long r = 12;   /* tg radius */
      /* get image */
      const mxArray* mx_im = prhs[0];
      matrix<> im = to_matrix(mx_im);
      /* get number of textons */
      unsigned long K = static_cast<unsigned long>(*mxGetPr(prhs[1]));
      /* get number of iterations */
      unsigned long n_iter = static_cast<unsigned long>(*mxGetPr(prhs[2]));
      /* get clustering problem size */
      double subsampling = *mxGetPr(prhs[3]);
      /* set # procs */
      if (nrhs > 4) {
         unsigned long n_proc = static_cast<unsigned long>(*mxGetPr(prhs[4]));
         thread::processors(n_proc);
      }  
      /* compute texton filter set */
      cout << "computing filter set... ";
cout.flush();
clock_t time = clock();
      auto_collection< matrix<>, array_list< matrix<> > > filters = 
         lib_image::texton_filters(n_ori, sigma);
cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* compute textons */
      cout << "computing textons...\n";
cout.flush();
      auto_collection< matrix<>, array_list< matrix<> > > textons;
      matrix<unsigned long> assign = lib_image::textons(
         im, *filters, textons, K, n_iter, subsampling
      );
      cout << "computing tg... ";
cout.flush();
time = clock();
      auto_collection< matrix<>, array_list< matrix<> > > gradients = 
         lib_image::hist_gradient_2D(assign, r, n_ori);
cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";

      cout << "smoothing tg... ";
cout.flush();
time = clock();
      auto_collection< matrix<>, array_list< matrix<> > > gradients_smoothed = 
         lib_image::smooth_grad_2D(*gradients, 17.0/3.0, 1);
cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      cout << "combining tg... ";
cout.flush();
time = clock();
      auto_ptr< matrix<> > tg;
      auto_ptr< matrix<> > tg_ori;
      lib_image::combine_hist_gradients(*gradients_smoothed, tg, tg_ori);
cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
/*
cout << "combining tg... ";
cout.flush();
time = clock();
      auto_ptr< matrix<> > tg;
      auto_ptr< matrix<> > tg_ori;
      lib_image::combine_hist_gradients(*gradients, tg, tg_ori);
cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
*/
      /* return texton assignments */
      if (nlhs > 0)
         plhs[0] = to_mxArray(matrix<>(assign));
      if (nlhs > 1)
         plhs[1] = to_mxArray(lib_image::nonmax_oriented_2D(*tg,*tg_ori));
      if (nlhs > 2)
         plhs[2] = to_mxArray(*tg);
      if (nlhs > 3)
         plhs[3] = to_mxArray(*tg_ori);
      if (nlhs > 4) {
         for (unsigned long n = 0; n < n_ori; n++)
            plhs[4+n] = to_mxArray((*gradients_smoothed)[n]);
      }
      if (nlhs > int(4+n_ori)) {
         for (unsigned long n = 0; n < n_ori; n++)
            plhs[4+n+n_ori] = to_mxArray((*gradients)[n]);
      }

      if (nlhs > int(4+2*n_ori)) {
         unsigned long support = 12;
         double sigma_y = (support > 0) ? (static_cast<double>(support)/3) : 1;
         double sigma_x = sigma_y;
         plhs[4+2*n_ori] = to_mxArray(lib_image::gaussian_2D(sigma_x, sigma_y, 0, 0, false, support, support));
         plhs[4+2*n_ori+1] = to_mxArray(lib_image::gaussian_2D(sigma_x, sigma_y, 0, 1, false, support, support));
         plhs[4+2*n_ori+2] = to_mxArray(lib_image::gaussian_2D(sigma_x, sigma_y, 0, 2, false, support, support));      
      }

   } catch (ex_index_out_of_bounds& e) {
      cout << "index: " << e.index() << "\n";
      cout << e << "\n";
   } catch (exception& e) {
      cout << e << "\n";
   }
}
