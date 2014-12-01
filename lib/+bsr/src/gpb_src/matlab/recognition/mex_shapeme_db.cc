/*
 * Database of geometric blur features.
 */
#include <iostream>

#include "collections/array_list.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "lang/array.hh"
#include "lang/exceptions/exception.hh"
#include "lang/pointers/auto_ptr.hh"
#include "lang/pointers/safe_ptr.hh"
#include "math/math.hh"
#include "math/matrices/matrix.hh"
#include "vision/features/geometric_blur.hh"
#include "vision/recognition/classifiers/unigram_classifier.hh"
#include "vision/recognition/classifiers/bigram_classifier.hh"
#include "vision/recognition/databases/category_db.hh"
#include "vision/recognition/models/category.hh"
#include "vision/recognition/models/exemplar.hh"
#include "vision/recognition/models/ids/category_id.hh"
#include "vision/recognition/models/ids/exemplar_id.hh"
#include "vision/recognition/models/ids/feature_id.hh"

#include "mlearning/clustering/clusterers/covering/matrix_clusterer.hh"

#include <cstring>
#include <mex.h>

using collections::array_list;
using collections::list;
using collections::pointers::auto_collection;
using lang::array;
using lang::exceptions::exception;
using lang::pointers::auto_ptr;
using lang::pointers::safe_ptr;
using math::matrices::matrix;
using vision::features::geometric_blur;
using vision::recognition::classifiers::unigram_classifier;
using vision::recognition::classifiers::bigram_classifier;
using vision::recognition::databases::category_db;
using vision::recognition::models::category;
using vision::recognition::models::exemplar;
using vision::recognition::models::ids::category_id;

using mlearning::clustering::clusterers::covering::matrix_clusterer;

/********************************** 
 * Matlab array conversion routines.
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
mxArray* mxCreateScalarDouble (double value) {
   mxArray* pa = mxCreateDoubleMatrix(1, 1, mxREAL);
   *mxGetPr(pa) = value;
   return pa;
}

/* 
 * Convert an mxArray to a double array.
 */
double** toDouble(const mxArray *a) {
   int mrows = mxGetM(a);
   int ncols = mxGetN(a);
   double *arr = mxGetPr(a);
   double *v = new double[(mrows * ncols)];
   double **data = new double*[mrows];
   int r, c;
   for (r = 0; r < mrows; r++) {
      data[r] = v;
      v += ncols;
   }
   /* copy the data */
   for (r = 0; r < mrows; r++) {
      for (c = 0; c < ncols; c++) {
         data[r][c] = arr[(c*mrows) + r];
      }
   }
   return data;
}

/*
 * Convert an mxArray to an array<double>.
 */
array<double> toArrayDouble(const mxArray* a) {
   int mrows = mxGetM(a);
   int ncols = mxGetN(a);
   double *arr = mxGetPr(a);
   unsigned long size = static_cast<unsigned long>(mrows*ncols);
   array<double> data(size);
   for (unsigned long n = 0; n < size; n++)
      data[n] = arr[n];
   return data;
}

/*
 * Convert an array<double> to an mxArray.
 */
mxArray* toMxArray(const array<double>& a) {
   mxArray *arr = mxCreateDoubleMatrix(
      1, static_cast<int>(a.size()), mxREAL
   );
   double *data = mxGetPr(arr);
   for (unsigned long n = 0; n < a.size(); n++)
      data[n] = a[n];
   return arr;
}

/*
 * Convert a list of array<double>'s to a multi-dimensional mxArray.
 */
mxArray* toMxArray(const list< array<double> >& lst) {
   int mrows = static_cast<int>(lst.size());
   int ncols = ((mrows > 0) ? (lst.head().size()) : 0);
   mxArray* arr = mxCreateDoubleMatrix(mrows, ncols, mxREAL);
   double* data = mxGetPr(arr);
   int r = 0;
   for (list< array<double> >::iterator_t i(lst); i.has_next(); ) {
      array<double>& a = i.next();
      for (int c = 0; c < ncols; c++) {
         data[(c*mrows) + r] = a[c];
      }
      r++;
   }
   return arr;
}

/*
 * Convert a double array back to an mxArray.
 */
mxArray* toMxArray(double **data, int mrows, int ncols) {
   mxArray *a = mxCreateDoubleMatrix(mrows, ncols, mxREAL);
   double *arr = mxGetPr(a);
   /* copy the data */
   int r, c;
   for (r = 0; r < mrows; r++) {
      for (c = 0; c < ncols; c++) {
         arr[(c*mrows) + r] = data[r][c];
      }
   }
   return a;
}

/*
 * Delete a double array.
 */
void deleteDouble(double **data) {
   delete [] (data[0]);
   delete [] data;
}

/*
 * Print an array.
 */
void print(const array<double> a) {
   for (unsigned long n = 0; n < a.size(); n++)
      std::cout << a[n] << " ";
}
 
/******************************************************************/

/*
 * Database.
 */
static array_list<category> c_list;   /* pending category list */

/*
 * Shapeme vocabulary.
 */
static array_list< matrix<> >                                    s_gb_list;   /* gb features to cluster */
static auto_collection< matrix<>, array_list< matrix<> > > s_centroids; /* shapeme centroids */

/*
 * Shapeme frequency counts.
 */
matrix<> count_C_s;  /* count_C_s(C,s) = # times shapeme s occurs in category C   */
matrix<> count_C;    /* count_C(C)     = total # shapeme occurances in category C */

/*
 * Utility function - get id of nearest cluster center.
 */
unsigned long find_cluster(
   const matrix<>& m, 
   const array_list< matrix<> >& centroids) 
{
   unsigned long min_id = 0;
   matrix<> diff = centroids[0] - m;
   double min_dist = dot(diff, diff);
   for (unsigned long n = 1; n < centroids.size(); n++) {
      matrix<> curr_diff = centroids[n] - m;
      double curr_dist = dot(curr_diff, curr_diff);
      if (curr_dist < min_dist) {
         min_id = n;
         min_dist = curr_dist;
      }
   }
   return min_id;
}

/*
 * Matlab interface - extract unigram features from an mxArray.
 * Input:
 *
 * mx_features - (feature vector length) x n_features
 *               matrix of geometric blur descriptors
 *
 * mx_x_pos    - vector of length n_features
 *               x-coordinate for each feature
 *
 * mx_y_pos    - vector of length n_features
 *               y-coordinate for each feature
 *
 * Output: exemplar containing features.
 */
auto_ptr<exemplar> matlab_extract_features(
   const mxArray* mx_features,
   const mxArray* mx_x_pos,
   const mxArray* mx_y_pos)
{
   /* extract features */
   array<double> feature_data = toArrayDouble(mx_features);
   array<double> x_pos = toArrayDouble(mx_x_pos);
   array<double> y_pos = toArrayDouble(mx_y_pos);
   unsigned long fvec_len   = static_cast<unsigned long>(mxGetM(mx_features));
   unsigned long n_features = static_cast<unsigned long>(mxGetN(mx_features));
   auto_collection<geometric_blur> features(new list<geometric_blur>());
   for (unsigned long n = 0; n < n_features; n++) {
      array<double> fvec = feature_data.subarray(
         n*fvec_len, 
         (n+1)*fvec_len - 1
      );
      matrix<> m_fvec(fvec.size(), 1);
      for (unsigned long i = 0; i < fvec.size(); i++)
         m_fvec[i] = fvec[i];
      auto_ptr<geometric_blur> gb(
         new geometric_blur(x_pos[n], y_pos[n], 0, 1, m_fvec)
      );
      features->add(*gb);
      gb.release();
   }
   return auto_ptr<exemplar>(new exemplar(features));
}

/*
 * Matlab interface - add exemplar to database.
 */
void matlab_add_exemplar_to_db(
   auto_ptr<exemplar> e,   /* exemplar */
   const mxArray* mx_c_id) /* category id */
{
   unsigned long n_categories = c_list.size();
   unsigned long c_id = static_cast<unsigned long>(*mxGetPr(mx_c_id));
   for (unsigned long n = n_categories; n < (c_id+1); n++) {
      auto_ptr<category> c(new category());
      c_list.add(*c);
      c.release();
   }
   c_list[c_id].add_exemplar(e);
}

/*
 * Matlab interface - add exemplar to vocabulary.
 */
void matlab_add_exemplar_to_vocab(auto_ptr<exemplar> e) {
   const array_list<geometric_blur>& gb_list = e->gb_features();
   for (unsigned long n = 0; n < gb_list.size(); n++) {
      const matrix<>& f = gb_list[n].descriptor();
      s_gb_list.add(*(new matrix<>(f)));
   }
} 

/*
 * Matlab interface - build database.
 */
void matlab_build_db(const mxArray* mx_r) {
   /* compute shapeme vocabulary by clustering */
std::cout << "clustering into shapemes...\n";
   double r = static_cast<double>(*mxGetPr(mx_r));
   matrix_clusterer<> cover_clust(r);
   cover_clust.cluster(s_gb_list, s_centroids);
   unsigned long K = s_centroids->size();
std::cout << K << " shapemes computed\n";
   /* compute shapeme frequency counts for each category */
std::cout << "computing category shapeme frequency counts...\n";
   unsigned long n_categories = c_list.size();
   count_C_s = matrix<>::ones(n_categories, K);
   count_C   = matrix<>::zeros(n_categories, 1);
   count_C.fill(static_cast<double>(K));
   for (unsigned long n = 0; n < n_categories; n++) {
std::cout << "\ncategory " << n << "\n";
      const category& c = c_list[n];
      const array_list<exemplar>& e_list = c.exemplars();
      for (unsigned long n_e = 0; n_e < e_list.size(); n_e++) {
         const exemplar& e = e_list[n_e];
         const array_list<geometric_blur>& gb_list = e.gb_features();
         for (unsigned long n_f = 0; n_f < gb_list.size(); n_f++) {
            const matrix<>& f = gb_list[n_f].descriptor();
            unsigned long s_id = find_cluster(f, *s_centroids);
            count_C_s(n,s_id)++;
            count_C[n]++;
         }
      }
   }
}

/*
 * Matlab interface - perform query.
 *
 * mx_features   - (feature vector length) x n_features
 *                 matrix of geometric blur descriptors
 *
 * mx_x_pos      - vector of length n_features
 *                 x-coordinate for each feature
 *
 * mx_y_pos      - vector of length n_features
 *                 y-coordinate for each feature
 * 
 * Return:
 *
 * mx_results    - vector of categories ids (one for each classifier run)
 */
void matlab_query_db(
   const mxArray* mx_features, 
   const mxArray* mx_x_pos,
   const mxArray* mx_y_pos,
   mxArray*& mx_results)
{
   /* create query exemplar */
   auto_ptr<exemplar> e = matlab_extract_features(
      mx_features, 
      mx_x_pos,
      mx_y_pos
   );
   /* compute shapemes in exemplar by vector quantization */
   const array_list<geometric_blur>& gb_list = e->gb_features();
   unsigned long n_features = gb_list.size();
   array<unsigned long> f_ids(n_features);
   for (unsigned long n = 0; n < n_features; n++) {
      const matrix<>& f = gb_list[n].descriptor();
      f_ids[n] = find_cluster(f, *s_centroids);
   }
   /* naive bayes classifier */
   unsigned long n_categories = count_C_s.size(0);
   array<double> category_probs(n_categories);
   /* compute log probability for each category */
   for (unsigned long n = 0; n < n_categories; n++) {
      double p = 0;
      for (unsigned long n_f = 0; n_f < n_features; n_f++)
         p += math::log(count_C_s(n, f_ids[n_f])/(count_C[n]));
      category_probs[n] = p;
   }
   /* find category with highest log probability */
   unsigned long c_id = 0;
   double c_prob = category_probs[0];
   for (unsigned long n = 1; n < n_categories; n++) {
      if (category_probs[n] > c_prob) {
         c_id = n;
         c_prob = category_probs[n];
      }
   }
   /* report results (category id) */
   array<double> results(1);
   results[0] = c_id;
   mx_results = toMxArray(results);
}

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
      if (strcmp(fun_name, "add_exemplar_to_db") == 0) {
         /* add an exemplar to the database */
         const mxArray* mx_features   = prhs[1];
         const mxArray* mx_x_pos      = prhs[2];
         const mxArray* mx_y_pos      = prhs[3];
         const mxArray* mx_c_id       = prhs[4];
         auto_ptr<exemplar> e = matlab_extract_features(
            mx_features, 
            mx_x_pos,
            mx_y_pos
         );
         matlab_add_exemplar_to_db(e, mx_c_id);
      } else if (strcmp(fun_name, "add_exemplar_to_vocab") == 0) {
         const mxArray* mx_features   = prhs[1];
         const mxArray* mx_x_pos      = prhs[2];
         const mxArray* mx_y_pos      = prhs[3];
         auto_ptr<exemplar> e = matlab_extract_features(
            mx_features, 
            mx_x_pos,
            mx_y_pos
         );
         matlab_add_exemplar_to_vocab(e);
      } else if (strcmp(fun_name, "build_db") == 0) {
         const mxArray* mx_r          = prhs[1];
         matlab_build_db(mx_r);
      } else if (strcmp(fun_name, "query_db") == 0) {
         /* run query */
         const mxArray* mx_features   = prhs[1];
         const mxArray* mx_x_pos      = prhs[2];
         const mxArray* mx_y_pos      = prhs[3];
         matlab_query_db( 
            mx_features, 
            mx_x_pos,
            mx_y_pos,
            plhs[0]
         );
      } else {
         std::cout << "Invalid command.\n";
      }
   } catch (exception& e) {
      std::cout << e << "\n";
   }  
   delete [] fun_name;
}
