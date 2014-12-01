/*
 * Feature identity.
 */
#include "lang/pointers/safe_ptr.hh"
#include "vision/recognition/models/exemplar.hh"
#include "vision/recognition/models/ids/exemplar_id.hh"
#include "vision/recognition/models/ids/feature_id.hh"

namespace vision {
namespace recognition {
namespace models {
namespace ids {
/*
 * Imports.
 */
using lang::pointers::safe_ptr;
using vision::recognition::models::exemplar;

/*
 * Constructor.
 */
feature_id::feature_id(exemplar& e, unsigned long id)
 : _exemplar(e),
   _id(id)
{ }

/*
 * Copy constructor.
 */
feature_id::feature_id(const feature_id& f)
 : _exemplar(f._exemplar),
   _id(f._id)
{ }

/*
 * Destructor.
 */
feature_id::~feature_id() {
   /* do nothing */
}

/*
 * Get the exemplar to which the feature belongs.
 */
exemplar& feature_id::parent_exemplar() const {
   return _exemplar;
}

/*
 * Get the feature number within the exemplar.
 */
unsigned long feature_id::id() const {
   return _id;
}

/*
 * Compare to another feature identity.
 */
int feature_id::compare_to(const feature_id& f) const {
   /* compare exemplars */
   safe_ptr<const exemplar_id>   e_id =   _exemplar.exemplar_identity();
   safe_ptr<const exemplar_id> f_e_id = f._exemplar.exemplar_identity();
   int e_ref_compare = 
      (&_exemplar < &(f._exemplar))
       ? (-1)
       : ((&_exemplar > &(f._exemplar)) ? 1 : 0);
   int e_compare = 
      (e_id.get() == NULL) 
       ? ((f_e_id.get() == NULL) ? e_ref_compare : -1)
       : ((f_e_id.get() == NULL) ? 1 : e_id->compare_to(*f_e_id));
   /* compare features */
   int f_compare = 
      (_id < f._id)
       ? (-1)
       : ((_id > f._id) ? 1 : 0);
   return (e_compare == 0) ? f_compare : e_compare;
}

} /* namespace ids */
} /* namespace models */
} /* namespace recognition */
} /* namespace vision */
