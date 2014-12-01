/*
 * Pb. 
 */
#include "collections/pointers/auto_collection.hh"
#include "collections/array_list.hh"
#include "io/streams/cout.hh"
#include "io/streams/iomanip.hh"
#include "io/streams/ios.hh"
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
using io::streams::ios;
using io::streams::iomanip::setiosflags;
using io::streams::iomanip::setw;
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
 * Suppress border responses.
 */
void suppress_border(matrix<>& m, unsigned long r) {
   unsigned long size_x = m.size(0);
   unsigned long size_y = m.size(1);
   for (unsigned long x = 0; x < size_x; x++) {
      for (unsigned long y = 0; y < size_y; y++) {
         if (((x < r) || ((x + r) >= size_x)) ||
             ((y < r) || ((y + r) >= size_y)))
            m(x,y) = 0;
      }
   }
}


/*
 * Matlab interface.
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
   try {
      /* get arguments */
      matrix<> pb = to_matrix(prhs[0]);
      matrix<> pb_ori = to_matrix(prhs[1]);
      double nonmax_ori_tol = *mxGetPr(prhs[2]);
      /* nonmax suppress pb */
      matrix<> pb_nmax = lib_image::nonmax_oriented_2D(
         pb, pb_ori, nonmax_ori_tol
      );
      /* return */
      if (nlhs > 0)
         plhs[0] = to_mxArray(pb_nmax);
   } catch (ex_index_out_of_bounds& e) {
      cout << "index: " << e.index() << "\n";
      cout << e << "\n";
   } catch (exception& e) {
      cout << e << "\n";
   }
}
