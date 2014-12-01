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
static auto_ptr<category_db> c_db;     /* category database */
static array_list<category>  c_list;   /* pending category list */

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
void matlab_add_exemplar(
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
 * Matlab interface - build database.
 */
void matlab_build_db() {
   /* clear the current database */
   c_db.reset(NULL);
   /* link all categories to the root category */
   auto_ptr<category> c_root(new category());
   while (!(c_list.is_empty())) {
      auto_ptr<category> c(&(c_list.remove_head()));
      c_root->add_subcategory(c);
   }
   /* build the new database */
   c_db.reset(new category_db(c_root));
}

/*
 * Matlab interface - perform query.
 * Report the bigrams found given:
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
 * mx_num_nn     - number of nearest neighbors to find in unigram search
 * 
 * mx_item_limit - max number of leaf nodes to visit in unigram search
 *
 * Return:
 *
 * mx_results    - vector of categories ids (one for each classifier run)
 */
void matlab_query_db(
   const mxArray* mx_features, 
   const mxArray* mx_x_pos,
   const mxArray* mx_y_pos,
   const mxArray* mx_num_nn, 
   const mxArray* mx_item_limit,
   mxArray*& mx_results)
{
   /* initialize classifier(s) */
   unigram_classifier uclass(   
      static_cast<unsigned long>(*mxGetPr(mx_num_nn)),
      static_cast<unsigned long>(*mxGetPr(mx_item_limit))
   );
/*
   bigram_classifier uclass(   
      static_cast<unsigned long>(*mxGetPr(mx_num_nn)),
      static_cast<unsigned long>(*mxGetPr(mx_item_limit))
   );
*/
   /* create query exemplar */
   auto_ptr<exemplar> e = matlab_extract_features(
      mx_features, 
      mx_x_pos,
      mx_y_pos
   );
   /* run query */
   safe_ptr<const category_id> c_id = uclass.classify(*c_db, *e);
   /* report results */
   array<double> results(1);
   if (c_id.get() != NULL)
      results[0] = static_cast<double>(c_id->id());
   else
      results[0] = 0;
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
      if (strcmp(fun_name, "add_exemplar") == 0) {
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
         matlab_add_exemplar(e, mx_c_id);
      } else if (strcmp(fun_name, "build_db") == 0) {
         /* build database */      
         matlab_build_db();
      } else if (strcmp(fun_name, "query_db") == 0) {
         /* check that database is non-null */
         if (c_db.get() == NULL) {
            std::cout << "cannot run query - database not yet built\n";
         } else {
            /* run query */
            const mxArray* mx_features   = prhs[1];
            const mxArray* mx_x_pos      = prhs[2];
            const mxArray* mx_y_pos      = prhs[3];
            const mxArray* mx_num_nn     = prhs[4];
            const mxArray* mx_item_limit = prhs[5];
            matlab_query_db( 
               mx_features, 
               mx_x_pos,
               mx_y_pos,
               mx_num_nn,
               mx_item_limit,
               plhs[0]
            );
         }
      } else {
         std::cout << "Invalid command.\n";
      }
   } catch (exception& e) {
      std::cout << e << "\n";
   }  
   delete [] fun_name;
}
