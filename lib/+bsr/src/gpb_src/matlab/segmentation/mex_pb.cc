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
      /* parameters */
      unsigned long n_ori       = 8;         /* number of orientations */
      unsigned long num_L_bins  = 25;        /* # bins for bg */
      unsigned long num_a_bins  = 25;        /* # bins for cg_a */
      unsigned long num_b_bins  = 25;        /* # bins for cg_b */
      double bg_smooth_sigma    = 0.1;       /* bg histogram smoothing sigma */
      double cg_smooth_sigma    = 0.05;      /* cg histogram smoothing sigma */
      unsigned long r_bg        = 6;         /* radius for bg */
      unsigned long r_cg        = 12;        /* radius for cg */
      unsigned long r_tg        = 12;        /* radius for tg */
      unsigned long border      = 15;        /* border pixels */
      double sigma_tg_filt      = 4.0;       /* sigma for tg filters */
      double sigma_tg_smooth    = 17.0/3.0;  /* sigma for tg smoothing */
      double epsilon            = 1;         /* epsilon for tg smoothing */
      double nonmax_ori_tol     = M_PIl/8;   /* tolerance for oriented nonmax */
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
      /* compute bg */
      cout << setw(40) << "computing bg ";
      cout.flush();
      time = clock();
      auto_collection< matrix<>, array_list< matrix<> > > bgs = 
         lib_image::hist_gradient_2D(Lq, r_bg, n_ori, bg_smooth_kernel);
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* compute cg */
      cout << setw(40) << "computing cg_a ";
      cout.flush();
      time = clock();
      auto_collection< matrix<>, array_list< matrix<> > > cgs_a = 
         lib_image::hist_gradient_2D(aq, r_cg, n_ori, cga_smooth_kernel);
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      cout << setw(40) << "computing cg_b ";
      cout.flush();
      time = clock();
      auto_collection< matrix<>, array_list< matrix<> > > cgs_b = 
         lib_image::hist_gradient_2D(bq, r_cg, n_ori, cgb_smooth_kernel);
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* add cg components */
      cout << setw(40) << "computing cg ";
      cout.flush();
      time = clock();
      auto_collection< matrix<>, array_list< matrix<> > > cgs(
         new array_list< matrix<> >()
      );
      for (unsigned long n = 0; n < n_ori; n++)
         cgs->add(*(new matrix<>(((*cgs_a)[n] + (*cgs_b)[n])/2.0)));
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* normalize bg, cg */
      cout << setw(40) << "normalizing bg, cg ";
      cout.flush();
      time = clock();
      double bg_max = 0;
      double cg_max = 0;
      for (unsigned long n = 0; n < n_ori; n++) {
         suppress_border((*bgs)[n], r_bg);
         suppress_border((*cgs)[n], r_cg);
      }
      for (unsigned long n = 0; n < n_ori; n++) {
         double bmax = max((*bgs)[n]);
         double cmax = max((*cgs)[n]);
         if (bmax > bg_max) { bg_max = bmax; }
         if (cmax > cg_max) { cg_max = cmax; }
      }
      for (unsigned long n = 0; n < n_ori; n++) {
         (*bgs)[n] /= bg_max;
         (*cgs)[n] /= cg_max;
      }
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* combine bg */
      cout << setw(40) << "combining bg ";
      cout.flush();
      time = clock();
      auto_ptr< matrix<> > bg;
      auto_ptr< matrix<> > bg_ori;
      lib_image::combine_hist_gradients(*bgs, bg, bg_ori);
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* combine cg */
      cout << setw(40) << "combining cg ";
      cout.flush();
      time = clock();
      auto_ptr< matrix<> > cg;
      auto_ptr< matrix<> > cg_ori;
      lib_image::combine_hist_gradients(*cgs, cg, cg_ori);
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* compute texton filter set */
      cout << setw(40) << "computing filter set for textons ";
      cout.flush();
      time = clock();
      auto_collection< matrix<>, array_list< matrix<> > > filters = 
         lib_image::texton_filters(n_ori, sigma_tg_filt);
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* compute textons */
      cout << setw(40) << "computing textons ";
      cout.flush();
      time = clock();
      auto_collection< matrix<>, array_list< matrix<> > > textons;
      matrix<unsigned long> t_assign = lib_image::textons(
         gray, *filters, textons 
      );
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* compute tg */
      cout << setw(40) << "computing tg ";
      cout.flush();
      time = clock();
      auto_collection< matrix<>, array_list< matrix<> > > tgs = 
         lib_image::hist_gradient_2D(t_assign, r_tg, n_ori);
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* normalize tg */
      cout << setw(40) << "normalizing tg ";
      cout.flush();
      time = clock();
      double tg_max = 0;
      for (unsigned long n = 0; n < n_ori; n++) {
         suppress_border((*tgs)[n], r_tg);
         double t_max = max((*tgs)[n]);
         if (t_max > tg_max) { tg_max = t_max; }
      }
      for (unsigned long n = 0; n < n_ori; n++)
         (*tgs)[n] /= tg_max; 
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* smooth tg */ 
      cout << setw(40) << "smoothing tg ";
      cout.flush();
      time = clock();
      auto_collection< matrix<>, array_list< matrix<> > > tgs_smoothed = 
         lib_image::smooth_grad_2D(*tgs, sigma_tg_smooth, epsilon);
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* normalize tg */
      cout << setw(40) << "normalizing smoothed tg ";
      cout.flush();
      time = clock();
      tg_max = 0;
      for (unsigned long n = 0; n < n_ori; n++) {
         double t_max = max((*tgs_smoothed)[n]);
         if (t_max > tg_max) { tg_max = t_max; }
      }
      for (unsigned long n = 0; n < n_ori; n++)
         (*tgs_smoothed)[n] /= tg_max; 
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* combine tg */
      cout << setw(40) << "combining tg ";
      cout.flush();
      time = clock();
      auto_ptr< matrix<> > tg;
      auto_ptr< matrix<> > tg_ori;
      lib_image::combine_hist_gradients(*tgs_smoothed, tg, tg_ori);
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* compute pb */
      cout << setw(40) << "computing pb ";
      cout.flush();
      time = clock();
      auto_collection< matrix<>, array_list< matrix<> > > pbs(
         new array_list< matrix<> >()
      );
      for (unsigned long n = 0; n < n_ori; n++) {
         matrix<> p = 
            ((*bgs)[n])/0.13*0.31 +
            ((*cgs)[n])/0.077*0.44 +
            ((*tgs_smoothed)[n])/0.063*0.063;
         p = ~(exp(-p + 3.08) + 1.0);
         unsigned long n_inds = p.size();
         for (unsigned long ind = 0; ind < n_inds; ind++) {
            if (p[ind] < 0) { p[ind] = 0; }
            if (p[ind] > 1) { p[ind] = 1; }
         }
         pbs->add(*(new matrix<>(p)));
      }
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* combine pb */
      cout << setw(40) << "combining pb ";
      cout.flush();
      time = clock();
      auto_ptr< matrix<> > pb;
      auto_ptr< matrix<> > pb_ori;
      lib_image::combine_hist_gradients(*pbs, pb, pb_ori);
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* nonmax suppress pb */
      cout << setw(40) << "nonmax suppressing pb ";
      cout.flush();
      time = clock();
      matrix<> pb_nmax = lib_image::nonmax_oriented_2D(
         *pb, *pb_ori, nonmax_ori_tol
      );
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";

      /* compute bg+cg+tg */
      cout << setw(40) << "computing bg+cg+tg ";
      cout.flush();
      time = clock();
      auto_collection< matrix<>, array_list< matrix<> > > pb_sums(
         new array_list< matrix<> >()
      );
      for (unsigned long n = 0; n < n_ori; n++) {
         matrix<> p = (((*bgs)[n]) + ((*cgs)[n]) + ((*tgs_smoothed)[n]))/3.0;
         pb_sums->add(*(new matrix<>(p)));
      }
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* combine bg+cg+tg */
      cout << setw(40) << "combining bg+cg+tg ";
      cout.flush();
      time = clock();
      auto_ptr< matrix<> > pb_sum;
      auto_ptr< matrix<> > pb_sum_ori;
      lib_image::combine_hist_gradients(*pb_sums, pb_sum, pb_sum_ori);
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      /* nonmax suppress bg+cg+tg */
      cout << setw(40) << "nonmax suppressing bg+cg+tg ";
      cout.flush();
      time = clock();
      matrix<> pb_sum_nmax = lib_image::nonmax_oriented_2D(
         *pb_sum, *pb_sum_ori, nonmax_ori_tol
      );
      cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";

      /* return */
      if (nlhs > 0)
         plhs[0] = to_mxArray(lib_image::border_trim_2D(*bg, border));
      if (nlhs > 1)
         plhs[1] = to_mxArray(lib_image::border_trim_2D(*cg, border));
      if (nlhs > 2)
         plhs[2] = to_mxArray(lib_image::border_trim_2D(*tg, border));
      if (nlhs > 3)
         plhs[3] = to_mxArray(lib_image::border_trim_2D(*pb_sum, border));
      if (nlhs > 4)
         plhs[4] = to_mxArray(lib_image::border_trim_2D(pb_sum_nmax, border));
      for (unsigned long n = 0; n < n_ori; n++) {
         if (nlhs > (static_cast<int>(n)+5))
            plhs[n+5] = to_mxArray(lib_image::border_trim_2D((*pb_sums)[n], border));
      }
   } catch (ex_index_out_of_bounds& e) {
      cout << "index: " << e.index() << "\n";
      cout << e << "\n";
   } catch (exception& e) {
      cout << e << "\n";
   }
}
