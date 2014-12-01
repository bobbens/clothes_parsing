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
      /* parameters - binning and smoothing */
      unsigned long n_ori       = 8;                     /* number of orientations */
      unsigned long num_L_bins  = 25;                    /* # bins for bg */
      unsigned long num_a_bins  = 25;                    /* # bins for cg_a */
      unsigned long num_b_bins  = 25;                    /* # bins for cg_b */
      double bg_smooth_sigma    = 0.1;                   /* bg histogram smoothing sigma */
      double cg_smooth_sigma    = 0.05;                  /* cg histogram smoothing sigma */
      unsigned long border      = 30;                    /* border pixels */
      double sigma_tg_filt_sm   = 2.0;                   /* sigma for small tg filters */
      double sigma_tg_filt_lg   = math::sqrt(2) * 2.0;   /* sigma for large tg filters */
      /* parameters - radii */
      unsigned long n_bg = 1;
      unsigned long n_cg = 1;
      unsigned long n_tg = 1;
      unsigned long r_bg[] = { 10 };
      unsigned long r_cg[] = { 20 };
      unsigned long r_tg[] = { 20 };
      /* compute bg histogram smoothing kernel */
      matrix<> bg_smooth_kernel =
         lib_image::gaussian(bg_smooth_sigma*num_L_bins);
      matrix<> cga_smooth_kernel =
         lib_image::gaussian(cg_smooth_sigma*num_a_bins);
      matrix<> cgb_smooth_kernel =
         lib_image::gaussian(cg_smooth_sigma*num_b_bins);
      /* get image */
      matrix<> L = to_matrix(prhs[0]);
      matrix<> a = to_matrix(prhs[1]);
      matrix<> b = to_matrix(prhs[2]);
      /* mirror border */
      L = lib_image::border_mirror_2D(L, border);
      a = lib_image::border_mirror_2D(a, border);
      b = lib_image::border_mirror_2D(b, border);
      /* convert to grayscale */
      clock_t start_time = clock();
      clock_t time = start_time;
      cout << setiosflags(ios::left);
      cout << setw(40) << "converting RGB to grayscale ";
      cout.flush();
      matrix<> gray = lib_image::grayscale(L,a,b);
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* gamma correct */
      lib_image::rgb_gamma_correct(L,a,b,2.5);
      /* convert to Lab */
      cout << setw(40) << "converting RGB to Lab ";
      cout.flush();
      time = clock();
      lib_image::rgb_to_lab(L,a,b);
      lib_image::lab_normalize(L,a,b);
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* quantize color channels  */
      cout << setw(40) << "quantizing color channels ";
      cout.flush();
      time = clock();
      matrix<unsigned long> Lq = lib_image::quantize_values(L, num_L_bins);
      matrix<unsigned long> aq = lib_image::quantize_values(a, num_a_bins);
      matrix<unsigned long> bq = lib_image::quantize_values(b, num_b_bins);
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* compute texton filter set */
      cout << setw(40) << "computing filter set for textons ";
      cout.flush();
      time = clock();
      auto_collection< matrix<>, array_list< matrix<> > > filters_small = 
         lib_image::texton_filters(n_ori, sigma_tg_filt_sm);
      auto_collection< matrix<>, array_list< matrix<> > > filters_large = 
         lib_image::texton_filters(n_ori, sigma_tg_filt_lg);
      array_list< matrix<> > filters;
      filters.add(*filters_small);
      filters.add(*filters_large);
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* compute textons */
      cout << setw(40) << "computing textons ";
      cout.flush();
      time = clock();
      auto_collection< matrix<>, array_list< matrix<> > > textons;
      matrix<unsigned long> t_assign = lib_image::textons(
         gray, filters, textons, 64 
      );
      t_assign = matrix<unsigned long>(
         lib_image::border_mirror_2D(
            lib_image::border_trim_2D(matrix<>(t_assign), border), border
         )
      );
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* return textons */
      plhs[0] = to_mxArray(lib_image::border_trim_2D(matrix<>(t_assign), border));
      unsigned long count = 1;
      /* compute bg at each radius */
      for (unsigned long rnum = 0; rnum < n_bg; rnum++) {
         /* compute bg */
         cout << setw(40) << "computing bg " << "r = " << r_bg[rnum] << " ";
         cout.flush();
         time = clock();
         auto_collection< matrix<>, array_list< matrix<> > > bgs = 
            lib_image::hist_gradient_2D(Lq, r_bg[rnum], n_ori, bg_smooth_kernel);
         cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
         /* return bg */
         mxArray* m = mxCreateCellMatrix(static_cast<int>(n_ori), 1);
         for (unsigned long n = 0; n < n_ori; n++)
            mxSetCell(m, static_cast<int>(n), to_mxArray(lib_image::border_trim_2D((*bgs)[n], border)));
         plhs[count++] = m;
      }
      /* compute cga at each radius */
      for (unsigned long rnum = 0; rnum < n_cg; rnum++) {
         /* compute cga */
         cout << setw(40) << "computing cg_a " << "r = " << r_cg[rnum] << " ";
         cout.flush();
         time = clock();
         auto_collection< matrix<>, array_list< matrix<> > > cgs_a = 
            lib_image::hist_gradient_2D(aq, r_cg[rnum], n_ori, cga_smooth_kernel);
         cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
         /* return cga */
         mxArray* m = mxCreateCellMatrix(static_cast<int>(n_ori), 1);
         for (unsigned long n = 0; n < n_ori; n++)
            mxSetCell(m, static_cast<int>(n), to_mxArray(lib_image::border_trim_2D((*cgs_a)[n], border)));
         plhs[count++] = m;
      }
      /* compute cgb at each radius */
      for (unsigned long rnum = 0; rnum < n_cg; rnum++) {
         /* compute cgb */
         cout << setw(40) << "computing cg_b " << "r = " << r_cg[rnum] << " ";
         cout.flush();
         time = clock();
         auto_collection< matrix<>, array_list< matrix<> > > cgs_b = 
            lib_image::hist_gradient_2D(bq, r_cg[rnum], n_ori, cgb_smooth_kernel);
         cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
         /* return cgb */
         mxArray* m = mxCreateCellMatrix(static_cast<int>(n_ori), 1);
         for (unsigned long n = 0; n < n_ori; n++)
            mxSetCell(m, static_cast<int>(n), to_mxArray(lib_image::border_trim_2D((*cgs_b)[n], border)));
         plhs[count++] = m;
      }
      /* compute tg at each radius */
      for (unsigned long rnum = 0; rnum < n_tg; rnum++) {
         /* compute tg */
         cout << setw(40) << "computing tg " << "r = " << r_tg[rnum] << " ";
         cout.flush();
         time = clock();
         auto_collection< matrix<>, array_list< matrix<> > > tgs = 
            lib_image::hist_gradient_2D(t_assign, r_tg[rnum], n_ori);
         cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
         /* return */
         mxArray* m = mxCreateCellMatrix(static_cast<int>(n_ori), 1);
         for (unsigned long n = 0; n < n_ori; n++)
            mxSetCell(m, static_cast<int>(n), to_mxArray(lib_image::border_trim_2D((*tgs)[n], border)));
         plhs[count++] = m;
      }
   } catch (ex_index_out_of_bounds& e) {
      cout << "index: " << e.index() << "\n";
      cout << e << "\n";
   } catch (exception& e) {
      cout << e << "\n";
   }
}
