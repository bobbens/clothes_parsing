/*
 * Database of geometric blur features.
 */
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "lang/array.hh"
#include "lang/exceptions/exception.hh"
#include "lang/exceptions/ex_index_out_of_bounds.hh"
#include "math/matrices/matrix.hh"
#include "mlearning/clustering/clusterers/covering/matrix_clusterer.hh"
#include "mlearning/clustering/clusterers/general/recursive_clusterer.hh"
#include "mlearning/clustering/clusterers/kmeans/matrix_clusterer.hh"
#include "mlearning/clustering/metrics/matrix_metrics.hh"

#include <iostream>
#include <cstring>
#include <mex.h>

using collections::list;
using collections::pointers::auto_collection;
using lang::array;
using lang::exceptions::exception;
using lang::exceptions::ex_index_out_of_bounds;
using math::matrices::matrix;
using mlearning::clustering::clusterers::general::recursive_clusterer;
using mlearning::clustering::metrics::matrix_metrics;

using namespace mlearning::clustering::clusterers;

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

/******************************************************************/

static auto_collection< matrix<>, list< matrix<> > > items(
   new list< matrix<> >()
);

/*
 * Matlab interface.
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
   /* get function to call */
   char* fun_name = mexGetString(prhs[0]);
   try {
      /* branch to appropriate function */
      if (strcmp(fun_name, "display") == 0) {
         const mxArray* a = prhs[1];
         matrix<> m = to_matrix(a);
         std::cout << m << "\n";
         plhs[0] = to_mxArray(m);
      } else if (strcmp(fun_name, "add_item") == 0) {
         const mxArray* a = prhs[1];
         matrix<> m = to_matrix(a);
         items->add(*(new matrix<>(m)));
      } else if (strcmp(fun_name, "clear") == 0) {
         items.reset(new list< matrix<> >());
      } else if (strcmp(fun_name, "kmeans_cluster_L2") == 0) {
         unsigned long K = static_cast<unsigned long>(*mxGetPr(prhs[1]));
         unsigned long M = static_cast<unsigned long>(*mxGetPr(prhs[2]));
         kmeans::matrix_clusterer<> km(K, M);
         array<unsigned long> assign = km.cluster(*items);
         unsigned long n_items = assign.size();
         matrix<> m(n_items, 1);
         for (unsigned long n = 0; n < n_items; n++)
            m[n] = static_cast<double>(assign[n]);
         plhs[0] = to_mxArray(m);
      }  else if (strcmp(fun_name, "kmeans_cluster_L1") == 0) {
         unsigned long K = static_cast<unsigned long>(*mxGetPr(prhs[1]));
         unsigned long M = static_cast<unsigned long>(*mxGetPr(prhs[2]));
         kmeans::matrix_clusterer<> km(K, M, matrix_metrics<>::L1_metric());
         array<unsigned long> assign = km.cluster(*items);
         unsigned long n_items = assign.size();
         matrix<> m(n_items, 1);
         for (unsigned long n = 0; n < n_items; n++)
            m[n] = static_cast<double>(assign[n]);
         plhs[0] = to_mxArray(m);
      } else if (strcmp(fun_name, "kmeans_cluster_L2_rec") == 0) {
         unsigned long K = static_cast<unsigned long>(*mxGetPr(prhs[1]));
         unsigned long M = static_cast<unsigned long>(*mxGetPr(prhs[2]));
         unsigned long B = static_cast<unsigned long>(*mxGetPr(prhs[3]));
         kmeans::matrix_clusterer<> km(K, M);
         recursive_clusterer< matrix<> > rc(km, B); 
         array<unsigned long> assign = rc.cluster(*items);
         unsigned long n_items = assign.size();
         matrix<> m(n_items, 1);
         for (unsigned long n = 0; n < n_items; n++)
            m[n] = static_cast<double>(assign[n]);
         plhs[0] = to_mxArray(m);
      } else if (strcmp(fun_name, "kmeans_cluster_L1_rec") == 0) {
         unsigned long K = static_cast<unsigned long>(*mxGetPr(prhs[1]));
         unsigned long M = static_cast<unsigned long>(*mxGetPr(prhs[2]));
         unsigned long B = static_cast<unsigned long>(*mxGetPr(prhs[3]));
         kmeans::matrix_clusterer<> km(K, M, matrix_metrics<>::L1_metric());
         recursive_clusterer< matrix<> > rc(km, B); 
         array<unsigned long> assign = rc.cluster(*items);
         unsigned long n_items = assign.size();
         matrix<> m(n_items, 1);
         for (unsigned long n = 0; n < n_items; n++)
            m[n] = static_cast<double>(assign[n]);
         plhs[0] = to_mxArray(m);
      } else if (strcmp(fun_name, "cover_cluster_L2") == 0) {
         double r = *mxGetPr(prhs[1]);
         covering::matrix_clusterer<> km(r);
         array<unsigned long> assign = km.cluster(*items);
         unsigned long n_items = assign.size();
         matrix<> m(n_items, 1);
         for (unsigned long n = 0; n < n_items; n++)
            m[n] = static_cast<double>(assign[n]);
         plhs[0] = to_mxArray(m);
      }  else if (strcmp(fun_name, "cover_cluster_L1") == 0) {
         double r = *mxGetPr(prhs[1]);
         covering::matrix_clusterer<> km(r, 0, matrix_metrics<>::L1_metric());
         array<unsigned long> assign = km.cluster(*items);
         unsigned long n_items = assign.size();
         matrix<> m(n_items, 1);
         for (unsigned long n = 0; n < n_items; n++)
            m[n] = static_cast<double>(assign[n]);
         plhs[0] = to_mxArray(m);
      } else if (strcmp(fun_name, "cover_cluster_L2_rec") == 0) {
         double r = *mxGetPr(prhs[1]);
         unsigned long B = static_cast<unsigned long>(*mxGetPr(prhs[2]));
         covering::matrix_clusterer<> km(r);
         recursive_clusterer< matrix<> > rc(km, B); 
         array<unsigned long> assign = rc.cluster(*items);
         unsigned long n_items = assign.size();
         matrix<> m(n_items, 1);
         for (unsigned long n = 0; n < n_items; n++)
            m[n] = static_cast<double>(assign[n]);
         plhs[0] = to_mxArray(m);
      } else if (strcmp(fun_name, "cover_cluster_L1_rec") == 0) {
         double r = *mxGetPr(prhs[1]);
         unsigned long B = static_cast<unsigned long>(*mxGetPr(prhs[2]));
         covering::matrix_clusterer<> km(r, 0, matrix_metrics<>::L1_metric());
         recursive_clusterer< matrix<> > rc(km, B); 
         array<unsigned long> assign = rc.cluster(*items);
         unsigned long n_items = assign.size();
         matrix<> m(n_items, 1);
         for (unsigned long n = 0; n < n_items; n++)
            m[n] = static_cast<double>(assign[n]);
         plhs[0] = to_mxArray(m);
      } else {
         std::cout << "Invalid command.\n";
      }
   } catch (ex_index_out_of_bounds& e) {
      std::cout << e << ": " << e.index() << "\n";
   } catch (exception& e) {
      std::cout << e << "\n";
   }  
   delete [] fun_name;
}
