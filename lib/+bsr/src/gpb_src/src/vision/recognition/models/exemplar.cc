/*
 * Exemplar.
 */
#include "collections/abstract/collection.hh"
#include "collections/array_list.hh"
#include "collections/pointers/auto_collection.hh"
#include "lang/iterators/iterator.hh"
#include "lang/null.hh"
#include "lang/pointers/auto_ptr.hh"
#include "lang/pointers/safe_ptr.hh"
#include "vision/features/geometric_blur.hh"
#include "vision/recognition/models/exemplar.hh"
#include "vision/recognition/models/ids/exemplar_id.hh"
#include "vision/recognition/models/ids/feature_id.hh"

namespace vision {
namespace recognition {
namespace models {
/*
 * Imports.
 */
using collections::abstract::collection;
using collections::array_list;
using collections::pointers::auto_collection;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;
using lang::pointers::safe_ptr;
using vision::features::geometric_blur;
using vision::recognition::models::ids::exemplar_id;
using vision::recognition::models::ids::feature_id;

/*
 * Constructor.
 * Create an exemplar from the given features.
 * Set the exemplar's identity.
 */
exemplar::exemplar(
   auto_collection<geometric_blur> gb_feats, 
   auto_ptr<exemplar_id> e_id)
 : _gb_features(NULL),
   _e_id(e_id)
{
   /* create array lists of features */
   auto_ptr< array_list<geometric_blur> > arr_gb_feats(new array_list<geometric_blur>());
   /* set all gb features to indicate exemplar as parent */
   unsigned long n = 0;
   for (auto_ptr< iterator<geometric_blur> > i = gb_feats->iter_create(); i->has_next(); n++) {
      geometric_blur& gb = i->next();
      auto_ptr<vision::features::feature_id> f_id(new feature_id(*this, n));
      gb.feature_identity(f_id);
      arr_gb_feats->add(gb);
   }
   /* release old ownership of features */
   auto_ptr< collection<geometric_blur> > c_gb_feats(gb_feats.release());
   /* transfer feature ownership to exemplar */
   _gb_features.reset(arr_gb_feats.release());
}

/*
 * Copy constructor.
 * Perform a deep copy of the exemplar.
 * Assign the copy the specified id.
 */
exemplar::exemplar(const exemplar& e, auto_ptr<exemplar_id> e_id)
 : _gb_features(new array_list<geometric_blur>),
   _e_id(e_id)
{
   /* copy and add gb features */
   unsigned long n = 0;
   for (auto_ptr< iterator<geometric_blur> > i = e._gb_features->iter_create(); i->has_next(); n++) {
      geometric_blur& gb = i->next();
      auto_ptr<vision::features::feature_id> f_id(new feature_id(*this, n));
      auto_ptr<geometric_blur> gb_copy(new geometric_blur(gb, f_id));
      _gb_features->add(*gb_copy);
      gb_copy.release();
   }
}

/*
 * Destructor.
 */
exemplar::~exemplar() {
   /* do nothing */
}

/*
 * Get geometric blur features in exemplar.
 */
const array_list<geometric_blur>& exemplar::gb_features() const {
   return *_gb_features;
}

/*
 * Get exemplar identity.
 */
safe_ptr<const exemplar_id> exemplar::exemplar_identity() const {
   return safe_ptr<const exemplar_id>(_e_id.get());
}

/*
 * Set exemplar identity.
 * Return the new identity.
 */
safe_ptr<const exemplar_id> exemplar::exemplar_identity(auto_ptr<exemplar_id> e_id) {
   _e_id = e_id;
   return safe_ptr<const exemplar_id>(_e_id.get());
}

} /* namespace models */
} /* namespace recognition */
} /* namespace vision */
