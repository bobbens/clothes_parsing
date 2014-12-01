/*
 * OE. 
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
      /* set # procs */
//      if (thread::processors() < 4)
//         thread::processors(4);
      /* set filter parameters */
      unsigned long n_ori = 8;
      double sigma = 5;
      /* get image */
      const mxArray* mx_im = prhs[0];
      matrix<> im = to_matrix(mx_im);
      /* get whether to use rectification */
      bool use_rectified = (nrhs > 1);
      /* compute filter sets */
      auto_collection< matrix<>, array_list< matrix<> > > fe;
      auto_collection< matrix<>, array_list< matrix<> > > fo;
      lib_image::oe_filters(n_ori, sigma, fe, fo);
      /* compute oe */
      if (!use_rectified) {
         /* oe */
         auto_collection< matrix<>, array_list< matrix<> > > oes_even;
         auto_collection< matrix<>, array_list< matrix<> > > oes_odd;
         lib_image::edge_oe(im, *fe, *fo, oes_even, oes_odd);
         /* combine oe */
         auto_ptr< matrix<> > edge;
         auto_ptr< matrix<> > edge_phase;
         auto_ptr< matrix<> > edge_ori;
         lib_image::combine_edge_oe(
            *oes_even, *oes_odd, edge, edge_phase, edge_ori
         );
         /* nonmax */
         matrix<> nmax = lib_image::nonmax_oriented_2D(*edge, *edge_ori, M_PI_4l/2);
         /* connected components */
         matrix<unsigned long> labels = lib_image::connected_components(nmax);
         /* return edges */
         if (nlhs > 0)
            plhs[0] = to_mxArray(nmax);
         if (nlhs > 1)
            plhs[1] = to_mxArray(*edge);
         if (nlhs > 2)
            plhs[2] = to_mxArray(*edge_phase);
         if (nlhs > 3)
            plhs[3] = to_mxArray(*edge_ori);
         if (nlhs > 4)
            plhs[4] = to_mxArray(matrix<>(labels));
      } else {
         /* rectified oe */
         auto_collection< matrix<>, array_list< matrix<> > > oes_neg;
         auto_collection< matrix<>, array_list< matrix<> > > oes_pos;
         lib_image::edge_oe_rectified(im, *fo, oes_neg, oes_pos);
         /* combine oe */
         auto_ptr< matrix<> > edge;
         auto_ptr< matrix<> > edge_polarity;
         auto_ptr< matrix<> > edge_ori;
         lib_image::combine_edge_oe_rectified(
            *oes_neg, *oes_pos, edge, edge_polarity, edge_ori
         );
         /* nonmax */
         matrix<> nmax = lib_image::nonmax_oriented_2D(*edge, *edge_ori);
         /* connected components */
         matrix<unsigned long> labels = lib_image::connected_components(nmax);
         /* return edges */
         if (nlhs > 0)
            plhs[0] = to_mxArray(nmax);
         if (nlhs > 1)
            plhs[1] = to_mxArray(*edge);
         if (nlhs > 2)
            plhs[2] = to_mxArray(*edge_polarity);
         if (nlhs > 3)
            plhs[3] = to_mxArray(*edge_ori);
         if (nlhs > 4)
            plhs[4] = to_mxArray(matrix<>(labels));
      }
   } catch (ex_index_out_of_bounds& e) {
      cout << "index: " << e.index() << "\n";
      cout << e << "\n";
   } catch (exception& e) {
      cout << e << "\n";
   }
}
