/*
 * Bigram classifier.
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
#include "vision/recognition/classifiers/bigram_classifier.hh"
#include "vision/recognition/models/exemplar.hh"
#include "vision/recognition/models/ids/category_id.hh"
#include "vision/recognition/models/ids/exemplar_id.hh"
#include "vision/recognition/models/ids/feature_id.hh"

#include <iostream> /*FIXME: remove this*/

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
bigram_classifier::bigram_classifier(
   unsigned long num_nn,
   unsigned long item_limit,
   double max_dist_change,
   double max_angle_change,
   double min_dist)
 : _num_nn(num_nn),
   _item_limit(item_limit),
   _max_dist_change(max_dist_change),
   _max_angle_change(max_angle_change),
   _min_dist(min_dist)
{ }

/*
 * Copy constructor.
 */
bigram_classifier::bigram_classifier(const bigram_classifier& c)
 : _num_nn(c._num_nn),
   _item_limit(c._item_limit),
   _max_dist_change(c._max_dist_change),
   _max_angle_change(c._max_angle_change),
   _min_dist(c._min_dist)
{ }

/*
 * Destructor.
 */
bigram_classifier::~bigram_classifier() {
   /* do nothing */
}

/*
 * Classify exemplar.
 */
safe_ptr<const category_id> bigram_classifier::classify(
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
   gb_search_opts.set_distance_limit(0.4);   /* FIXME: should not be hardcoded */
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
   /* create linked arrays of matches */
   array_list<const feature_id> gb_train_feat_ids;
   array_list<geometric_blur> gb_train_feats;
   array_list<geometric_blur> gb_test_feats;
   unsigned long n = 0;
   for (array_list< list<geometric_blur> >::iterator_t i(*gb_nn_lists); i.has_next(); n++) {
      geometric_blur&       gb_test = gb_feats[n];
      list<geometric_blur>& nn_list = i.next();
      for (list<geometric_blur>::iterator_t i_nn(nn_list); i_nn.has_next(); ) {
         geometric_blur& gb_train = i_nn.next();
         const feature_id&  f_id = dynamic_typecast<const feature_id&>(*(gb_train.feature_identity()));
         gb_train_feat_ids.add(f_id);
         gb_train_feats.add(gb_train);
         gb_test_feats.add(gb_test);
      }
   }
   /* sort by training feature id */
   array_list<const feature_id> gb_train_feat_ids_sorted(gb_train_feat_ids);
   array_list<geometric_blur> gb_train_feats_sorted;
   array_list<geometric_blur> gb_test_feats_sorted;
   array<unsigned long> inds = gb_train_feat_ids_sorted.sort_idx();
   gb_train_feats.subarray(inds, gb_train_feats_sorted);
   gb_test_feats.subarray(inds, gb_test_feats_sorted);
   /* create neighbor lists for training features (invert feature mapping above) */
   array_list<geometric_blur> gb_unique_train_feats;
   auto_collection< list<geometric_blur>, array_list< list<geometric_blur> > > 
      inv_gb_nn_lists(new array_list< list<geometric_blur> >);
   gb_unique_train_feats.add(gb_train_feats_sorted[0]);
   auto_ptr< list<geometric_blur> > temp_list(new list<geometric_blur>());
   temp_list->add(gb_test_feats_sorted[0]);
   inv_gb_nn_lists->add(*temp_list);
   temp_list.release();
   for (n = 1; n < gb_train_feats_sorted.size(); n++) {
      unsigned long n_unique = gb_unique_train_feats.size();
      const feature_id& gb_prev_id = gb_train_feat_ids_sorted[n-1];
      const feature_id& gb_curr_id = gb_train_feat_ids_sorted[n];
      if (gb_prev_id.compare_to(gb_curr_id) == 0) {
         /* curr feature is same as prev */
         (*inv_gb_nn_lists)[n_unique-1].add(gb_test_feats_sorted[n]);
      } else {
         /* curr feature is diff from prev */
         gb_unique_train_feats.add(gb_train_feats_sorted[n]);
         auto_ptr< list<geometric_blur> > temp_gb_list(new list<geometric_blur>());
         temp_gb_list->add(gb_test_feats_sorted[n]);
         inv_gb_nn_lists->add(*temp_gb_list);
         temp_gb_list.release();
      }
   }
   /* compute number of unique start/stop inds */
   unsigned long n_unique_start_stop = 1;
   for (n = 1; n < gb_unique_train_feats.size(); n++) {
      geometric_blur& gb_prev = gb_unique_train_feats[n-1];
      geometric_blur& gb_curr = gb_unique_train_feats[n];
      const feature_id& f_id_prev = dynamic_typecast<const feature_id&>(*(gb_prev.feature_identity()));
      const feature_id& f_id_curr = dynamic_typecast<const feature_id&>(*(gb_curr.feature_identity()));
      const exemplar_id& e_id_prev = *(f_id_prev.parent_exemplar().exemplar_identity());
      const exemplar_id& e_id_curr = *(f_id_curr.parent_exemplar().exemplar_identity());
      if (e_id_prev.compare_to(e_id_curr) != 0)
         n_unique_start_stop++;
   }
   /* record exemplar start/stop inds */
   array<unsigned long> e_start_inds(n_unique_start_stop);
   array<unsigned long> e_end_inds(n_unique_start_stop);
   e_start_inds[0] = 0;
   unsigned long count = 1;
   for (n = 1; n < gb_unique_train_feats.size(); n++) {
      geometric_blur& gb_prev = gb_unique_train_feats[n-1];
      geometric_blur& gb_curr = gb_unique_train_feats[n];
      const feature_id& f_id_prev = dynamic_typecast<const feature_id&>(*(gb_prev.feature_identity()));
      const feature_id& f_id_curr = dynamic_typecast<const feature_id&>(*(gb_curr.feature_identity()));
      const exemplar_id& e_id_prev = *(f_id_prev.parent_exemplar().exemplar_identity());
      const exemplar_id& e_id_curr = *(f_id_curr.parent_exemplar().exemplar_identity());
      if (e_id_prev.compare_to(e_id_curr) != 0) {
         e_start_inds[count] = n;
         e_end_inds[count-1] = n - 1;
         count++;
      }
   }
   e_end_inds[n_unique_start_stop-1] = gb_unique_train_feats.size() - 1;
   /* initialize quality measure for each match */
   auto_collection< array<double>, array_list< array<double> > >
      inv_gb_match_quality(new array_list< array<double> >);
   for (array_list< list<geometric_blur> >::iterator_t i(*inv_gb_nn_lists); i.has_next(); ) {
      list<geometric_blur>& nn_list = i.next();
      inv_gb_match_quality->add(*(new array<double>(nn_list.size())));
   }
   /* compute match quality measure */
   double mean_quality = 0;
   unsigned long total_matches = 0;
   auto_collection< double, array_list<double> > qualities(new array_list<double>);
   for (unsigned long e_num = 0; e_num < e_start_inds.size(); e_num++) {
      /* get exemplar feature start/end indices */
      unsigned long e_start = e_start_inds[e_num];
      unsigned long e_end   = e_end_inds[e_num];
      /* check bigrams between features in exemplar */
      for (unsigned long f1_train_ind = e_start; f1_train_ind <= e_end; f1_train_ind++) {
         geometric_blur&       f1_train           = gb_unique_train_feats[f1_train_ind];
         list<geometric_blur>& f1_matches         = (*inv_gb_nn_lists)[f1_train_ind];
         array<double>&        f1_matches_quality = (*inv_gb_match_quality)[f1_train_ind];
         unsigned long f1_test_ind = 0;
         for (list<geometric_blur>::iterator_t i_f1(f1_matches); i_f1.has_next(); f1_test_ind++) {
            geometric_blur& f1_test = i_f1.next();
            /* compute quality of match based on strength of other bigrams that verify it */
            double quality = 0;
            for (unsigned long f2_train_ind = e_start; f2_train_ind <= e_end; f2_train_ind++) {
               geometric_blur&       f2_train   = gb_unique_train_feats[f2_train_ind];
               list<geometric_blur>& f2_matches = (*inv_gb_nn_lists)[f2_train_ind];
               if (f1_train_ind != f2_train_ind) {
                  for (list<geometric_blur>::iterator_t i_f2(f2_matches); i_f2.has_next(); ) {
                     geometric_blur& f2_test = i_f2.next();
                     if (&f1_test != &f2_test) {
                        /* check bigram match */
                        double q = this->verify_match(
                           f1_train, f2_train, f1_test, f2_test
                        );
                        quality += q;
                     }
                  }
               }
            }
            /* update match quality */
            f1_matches_quality[f1_test_ind] = quality;
            /* update mean quality */
            if (quality > 0) {
               double* tempdouble = new double;
               *tempdouble = quality;
               qualities->add(*tempdouble);
               mean_quality += quality;
               total_matches++;
            }
         }
      }
   }
   /* sort qualities */
   qualities->sort();
   /* display mean quality */
   mean_quality /= double(total_matches);
   std::cout << "mean quality     = " << mean_quality << "\n";
   //double required_quality = (*qualities)[9*(qualities->size()/10)];
   double required_quality = mean_quality;
   std::cout << "required quality = " << required_quality << "\n";
   /* initialize category id votes */
   unsigned long n_categories = max_id + 1;
   array<double> votes(n_categories);
   array< safe_ptr<const category_id> > ids(n_categories);
   /* vote bigram pairs for which both matches have high quality */
   for (unsigned long e_num = 0; e_num < e_start_inds.size(); e_num++) {
      /* get exemplar feature start/end indices */
      unsigned long e_start = e_start_inds[e_num];
      unsigned long e_end   = e_end_inds[e_num];
      /* get class id of exemplar */
      geometric_blur& gb = gb_unique_train_feats[e_start];
      const feature_id&  f_id = dynamic_typecast<const feature_id&>(*(gb.feature_identity()));
      const exemplar_id& e_id = *(f_id.parent_exemplar().exemplar_identity());
      const category_id& c_id = *(e_id.parent_category().category_identity());
      unsigned long id = c_id.id();
      /* check bigrams between features in exemplar */
      for (unsigned long f1_train_ind = e_start; f1_train_ind <= e_end; f1_train_ind++) {
         geometric_blur&       f1_train           = gb_unique_train_feats[f1_train_ind];
         list<geometric_blur>& f1_matches         = (*inv_gb_nn_lists)[f1_train_ind];
         array<double>&        f1_matches_quality = (*inv_gb_match_quality)[f1_train_ind];
         unsigned long f1_test_ind = 0;
         for (list<geometric_blur>::iterator_t i_f1(f1_matches); i_f1.has_next(); f1_test_ind++) {
            geometric_blur& f1_test    = i_f1.next();
            double          f1_quality = f1_matches_quality[f1_test_ind];
            if (f1_quality > required_quality) {
               for (unsigned long f2_train_ind = e_start; f2_train_ind <= e_end; f2_train_ind++) {
                  geometric_blur&       f2_train           = gb_unique_train_feats[f2_train_ind];
                  list<geometric_blur>& f2_matches         = (*inv_gb_nn_lists)[f2_train_ind];
                  array<double>&        f2_matches_quality = (*inv_gb_match_quality)[f2_train_ind];
                  unsigned long f2_test_ind = 0;
                  if (f1_train_ind != f2_train_ind) {
                     for (list<geometric_blur>::iterator_t i_f2(f2_matches); i_f2.has_next(); f2_test_ind++) {
                        geometric_blur& f2_test    = i_f2.next();
                        double          f2_quality = f2_matches_quality[f2_test_ind];
                        if ((f2_quality > required_quality) && (&f1_test != &f2_test)) {
                           /* check bigram match */
                           double q = this->verify_match(
                              f1_train, f2_train, f1_test, f2_test
                           );
                           if (q > 0) {
                              /* cast vote */
                              votes[id] += 1;
                              ids[id].reset(&c_id);
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   /* find category with the most votes */
   unsigned long best_category = 0;
   double best_score = votes[0];
   for (n = 0; n < n_categories; n++) {
      if (votes[n] > best_score) {
         best_category = n;
         best_score    = votes[n];
      }
std::cout << votes[n] << " ";
   }
std::cout << "\n";
   /* return id of best category */
   return ids[best_category];
}

double std_angle(double theta) {
   const double PI = 3.141592653589793115997963468544185161590576171875;
   while (theta < 0) 
      theta += 2*PI;
   while (theta >= 2*PI)
      theta -= 2*PI;
   return theta;
}
   
/*
 * Verify match.
 * Return match quality.
 */
double bigram_classifier::verify_match(
   const geometric_blur& f1_train, 
   const geometric_blur& f2_train, 
   const geometric_blur& f1_test, 
   const geometric_blur& f2_test)
{
   /* extract feature locations */
   double f1_train_x = f1_train.x();
   double f1_train_y = f1_train.y();
   double f2_train_x = f2_train.x();
   double f2_train_y = f2_train.y();
   double f1_test_x = f1_test.x();
   double f1_test_y = f1_test.y();
   double f2_test_x = f2_test.x();
   double f2_test_y = f2_test.y();
   /* compute distances */
   double dist_train_x = f1_train_x - f2_train_x;
   double dist_train_y = f1_train_y - f2_train_y;
   double dist_test_x = f1_test_x - f2_test_x;
   double dist_test_y = f1_test_y - f2_test_y;
   double dist_train = math::sqrt(dist_train_x*dist_train_x + dist_train_y*dist_train_y);
   double dist_test  = math::sqrt(dist_test_x*dist_test_x   + dist_test_y*dist_test_y);
   double dist_change = 
      (dist_train < dist_test)
       ? ((dist_test - dist_train)/dist_train)
       : ((dist_train - dist_test)/dist_test);
   /* compute angle change */
   double ori_train = math::atan2(dist_train_y, dist_train_x);
   double ori_test  = math::atan2(dist_test_y,  dist_test_x);
   double angle_changeA = std_angle(ori_train - ori_test);
   double angle_changeB = std_angle(ori_test  - ori_train);
   double angle_change = (angle_changeA < angle_changeB) ? angle_changeA : angle_changeB;
   /* compute match quality */
//   double match_quality = (2.0 - f1_train.distance_to(f1_test)) + (2.0 - f2_train.distance_to(f2_test));
   /* determine verification */
   bool is_verified = ((dist_change < _max_dist_change) &&
                       (angle_change < _max_angle_change) &&
                       (dist_train > _min_dist));
//   return (is_verified) ? match_quality : 0;
   return (is_verified) ? 1 : 0;
}

} /* namespace classifiers */
} /* namespace recognition */
} /* namespace vision */
