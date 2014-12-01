/*
 * Geometric Blur.
 */
#include "collections/abstract/collection.hh"
#include "collections/list.hh"
#include "interfaces/kd_treeable.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/matrices/functors/matrix_distance_functors.hh"
#include "math/matrices/functors/matrix_equal_functors.hh"
#include "math/matrices/functors/matrix_key_functors.hh"
#include "math/matrices/matrix.hh"
#include "vision/features/feature_descriptor.hh"
#include "vision/features/feature_id.hh"
#include "vision/features/geometric_blur.hh"
#include "vision/features/parameters/geometric_blur_params.hh"

namespace vision {
namespace features {
/*
 * Imports.
 */
using collections::abstract::collection;
using collections::list;
using interfaces::kd_tree_key;
using lang::exceptions::ex_invalid_argument;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;
using math::matrices::functors::matrix_distance_functors;
using math::matrices::functors::matrix_equal_functors;
using math::matrices::functors::matrix_key_functors;
using math::matrices::functors::matrix_L2_distance;
using math::matrices::functors::matrix_equal;
using math::matrices::functors::matrix_key_compare;
using math::matrices::functors::matrix_key_split;
using math::matrices::functors::matrix_key_distance;
using math::matrices::matrix;
using vision::features::parameters::geometric_blur_params;

/*
 * Constructor.
 * Create a geometric blur feature as specified.
 * FIXME: this should be removed in the future 
 * (once c++ implementation of computing features is done).
 */
geometric_blur::geometric_blur(
   double x_coord, 
   double y_coord, 
   double ori, 
   double scl,
   const matrix<>& desc,
   auto_ptr<feature_id> f_id,
   const geometric_blur_params& params)
 : feature_descriptor(
      x_coord, 
      y_coord, 
      ori, 
      scl,
      desc,
      f_id),
   _params(params)
{ }

/*
 * Copy constructor.
 * Assign the copy the specified id.
 */
geometric_blur::geometric_blur(
   const geometric_blur& gb,
   auto_ptr<feature_id> f_id)
 : feature_descriptor(gb, f_id),
   _params(gb._params)
{ }

/*
 * Destructor.
 */
geometric_blur::~geometric_blur() {
   /* do nothing */
}

/*
 * Get the parameters associated with the feature.
 */
const geometric_blur_params& geometric_blur::parameters() const {
   return _params;
}

/*
 * Distance from feature vector to half-space.
 */
double geometric_blur::distance_to(const kd_tree_key<>& k) const {
   const matrix_key_distance<>& f = matrix_key_functors<>::f_key_distance();
   return f(_descriptor, k);
}

/*
 * Distance between feature vectors.
 * Throw an invalid argument exception if the features are not comparable
 * (do not share the same parameters).
 */
double geometric_blur::distance_to(const geometric_blur& gb) const {
   /* check parameter compatibility */
   if (&_params != &(gb._params)) {
      throw ex_invalid_argument(
         "attempted comparison of incompatible geometric blur features"
      );
   }
   /* compute distance */
   const matrix_L2_distance<>& f = matrix_distance_functors<>::L2_distance();
   return f(_descriptor, gb._descriptor);
}

/*
 * Determine side of half-space in which feature vector lies.
 */
int geometric_blur::compare_to(const kd_tree_key<>& k) const {
   const matrix_key_compare<>& f = matrix_key_functors<>::f_key_compare();
   return f(_descriptor, k);
}

/*
 * Check equality of feature vectors.
 * Throw an invalid argument exception if the features are not comparable
 * (do not share the same parameters).
 */
bool geometric_blur::is_equal_to(const geometric_blur& gb) const {
   /* check parameter compatibility */
   if (&_params != &(gb._params)) {
      throw ex_invalid_argument(
         "attempted comparison of incompatible geometric blur features"
      );
   }
   /* perform comparison */
   const matrix_equal<>& f = matrix_equal_functors<>::f_equal();
   return f(_descriptor, gb._descriptor);
}

/*
 * Return dimensionality of geometric blur feature vectors in collection.
 * Throw an invalid argument exception if the items in the collection do 
 * not share the same parameters (and hence have the same dimensionality).
 * An invalid argument exception is also thrown if the collection is empty.
 */
unsigned long geometric_blur::dimensionality(const collection<geometric_blur>& c) {
   /* check that all features share the same parameters */
   geometric_blur::assert_nonempty_compatible(c);
   /* get dimensionality */
   auto_ptr< iterator<geometric_blur> > i = c.iter_create();
   const geometric_blur_params& p = i->next()._params;
   return p.descriptor_size();
}

/*
 * Compute key for splitting a collection of geometric blur features 
 * along the dimension of the feature vector with highest variance.  
 * Throw an invalid argument exception if the features do not share 
 * the same parameters (or the collection is empty).
 */
kd_tree_key<> geometric_blur::key_split(const collection<geometric_blur>& c) {
   /* check that all features share the same parameters */
   geometric_blur::assert_nonempty_compatible(c);
   /* build collection of feature vectors */
   list< matrix<> > descriptors;
   for (auto_ptr< iterator<geometric_blur> > i = c.iter_create(); i->has_next(); )
      descriptors.add(i->next()._descriptor);
   /* compute key */
   const matrix_key_split<>& f = matrix_key_functors<>::f_key_split();
   return f(descriptors);
}

/*
 * Assert that the given collection of geometric blur features is nonempty 
 * and that all features in the collection share the same parameters.
 * Throw an invalid argument exception if this assertation fails.
 */
void geometric_blur::assert_nonempty_compatible(const collection<geometric_blur>& c) {
   auto_ptr< iterator<geometric_blur> > i = c.iter_create();
   if (i->has_next()) {
      const geometric_blur_params& p = i->next()._params;
      while (i->has_next()) {
         if (&p != &(i->next()._params)) { 
            throw ex_invalid_argument(
               "not all geometric blur features in collection share the same parameters"
            );
         }
      }
   } else {
      throw ex_invalid_argument(
         "collection of geometric blur features is empty"
      );
   }
}

} /* namespace features */
} /* namespace vision */
