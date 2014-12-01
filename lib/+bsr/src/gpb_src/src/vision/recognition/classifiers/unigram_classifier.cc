/*
 * Unigram classifier.
 */
#include "collections/array_list.hh"
#include "collections/kd_tree.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "lang/array.hh"
#include "lang/pointers/auto_ptr.hh"
#include "lang/pointers/safe_ptr.hh"
#include "lang/typecasts/dynamic_typecast.hh"
#include "math/math.hh"
#include "vision/features/geometric_blur.hh"
#include "vision/recognition/databases/category_db.hh"
#include "vision/recognition/classifiers/unigram_classifier.hh"
#include "vision/recognition/models/exemplar.hh"
#include "vision/recognition/models/ids/category_id.hh"
#include "vision/recognition/models/ids/exemplar_id.hh"
#include "vision/recognition/models/ids/feature_id.hh"

namespace vision {
namespace recognition {
namespace classifiers {
/*
 * Imports.
 */
using collections::array_list;
using collections::kd_tree;
using collections::kd_tree_search_options;
using collections::list;
using collections::pointers::auto_collection;
using lang::array;
using lang::pointers::auto_ptr;
using lang::pointers::safe_ptr;
using lang::typecasts::dynamic_typecast;
using vision::features::geometric_blur;
using vision::recognition::databases::category_db;
using vision::recognition::models::exemplar;
using vision::recognition::models::ids::category_id;
using vision::recognition::models::ids::exemplar_id;
using vision::recognition::models::ids::feature_id;

/*
 * Constructor.
 * Specify parameters for the classifier.
 */
unigram_classifier::unigram_classifier(
   unsigned long num_nn, 
   unsigned long item_limit)
 : _num_nn(num_nn),
   _item_limit(item_limit)
{ }

/*
 * Copy constructor.
 */
unigram_classifier::unigram_classifier(const unigram_classifier& c)
 : _num_nn(c._num_nn),
   _item_limit(c._item_limit)
{ }

/*
 * Destructor.
 */
unigram_classifier::~unigram_classifier() {
   /* do nothing */
}

/*
 * Classify exemplar.
 */
safe_ptr<const category_id> unigram_classifier::classify(
   const category_db& c_db, 
   const exemplar& e)
{
   /* get list of features in exemplar */
   const array_list<geometric_blur>& gb_feats = e.gb_features();
   /* get feature search structure */
   const kd_tree<geometric_blur>& gb_db = c_db.gb_features();
   /* initialize parameters for finding nearest neighbors */
   kd_tree_search_options<> gb_search_opts(_num_nn);
   gb_search_opts.set_item_limit(_item_limit);
   /* find approximate nearest neighbors of features in database */
   auto_collection< list<geometric_blur>, array_list< list<geometric_blur> > > 
      gb_nn_lists(new array_list< list<geometric_blur> >);
   for (array_list<geometric_blur>::iterator_t i(gb_feats); i.has_next(); ) {
      /* grab the next feature */
      geometric_blur& gb = i.next();
      /* find nearest neighbors */
      list<geometric_blur> nn_list = gb_db.find_nns(gb, gb_search_opts);
      /* save nearest neighbor list */
      auto_ptr< list<geometric_blur> > nn_list_copy(new list<geometric_blur>(nn_list));
      gb_nn_lists->add(*nn_list_copy);
      nn_list_copy.release();
   }
   /* record maximum category id */
   unsigned long max_id = 0;
   for (array_list< list<geometric_blur> >::iterator_t i(*gb_nn_lists); i.has_next(); ) {
      list<geometric_blur>& nn_list = i.next();
      for (list<geometric_blur>::iterator_t i_nn(nn_list); i_nn.has_next(); ) {
         geometric_blur& gb = i_nn.next();
         const feature_id&  f_id = dynamic_typecast<const feature_id&>(*(gb.feature_identity()));
         const exemplar_id& e_id = *(f_id.parent_exemplar().exemplar_identity());
         const category_id& c_id = *(e_id.parent_category().category_identity());
         unsigned long id = c_id.id();
         if (id > max_id)
            max_id = id;
      }
   }
   /* initialize category id votes */
   unsigned long n_categories = max_id + 1;
   array<double> votes(n_categories);
   array< safe_ptr<const category_id> > ids(n_categories);
   /* record category id votes */
   for (array_list< list<geometric_blur> >::iterator_t i(*gb_nn_lists); i.has_next(); ) {
      list<geometric_blur>& nn_list = i.next();
      for (list<geometric_blur>::iterator_t i_nn(nn_list); i_nn.has_next(); ) {
         geometric_blur& gb = i_nn.next();
         const feature_id&  f_id = dynamic_typecast<const feature_id&>(*(gb.feature_identity()));
         const exemplar_id& e_id = *(f_id.parent_exemplar().exemplar_identity());
         const category_id& c_id = *(e_id.parent_category().category_identity());
         unsigned long id = c_id.id();
         /* weigh vote by number of features in exemplar */
         unsigned long n_e_features = f_id.parent_exemplar().gb_features().size();
         votes[id] += math::sqrt((double(1))/(double(n_e_features)));
         ids[id].reset(&c_id);
      }
   }
   /* find category with the most votes */
   unsigned long best_category = 0;
   double best_score = votes[0];
   for (unsigned long n = 1; n < n_categories; n++) {
      if (votes[n] > best_score) {
         best_category = n;
         best_score    = votes[n];
      }
   }
   /* return id of best category */
   return ids[best_category];
}

} /* namespace classifiers */
} /* namespace recognition */
} /* namespace vision */
