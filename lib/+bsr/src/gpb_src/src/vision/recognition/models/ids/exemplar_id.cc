/*
 * Exemplar identity.
 */
#include "lang/pointers/safe_ptr.hh"
#include "vision/recognition/models/category.hh"
#include "vision/recognition/models/ids/category_id.hh"
#include "vision/recognition/models/ids/exemplar_id.hh"

namespace vision {
namespace recognition {
namespace models {
namespace ids {
/*
 * Imports.
 */
using lang::pointers::safe_ptr;
using vision::recognition::models::category;

/*
 * Constructor.
 */
exemplar_id::exemplar_id(category& c, unsigned long id)
 : _category(c),
   _id(id)
{ }

/*
 * Copy constructor.
 */
exemplar_id::exemplar_id(const exemplar_id& e)
 : _category(e._category),
   _id(e._id)
{ }

/*
 * Destructor.
 */
exemplar_id::~exemplar_id() {
   /* do nothing */
}

/*
 * Get the category to which the exemplar belongs.
 */
category& exemplar_id::parent_category() const {
   return _category;
}

/*
 * Get the exemplar number within the category.
 */
unsigned long exemplar_id::id() const {
   return _id;
}

/*
 * Compare to another exemplar identity.
 */
int exemplar_id::compare_to(const exemplar_id& e) const {
   /* compare categories */
   safe_ptr<const category_id>   c_id =   _category.category_identity();
   safe_ptr<const category_id> e_c_id = e._category.category_identity();
   int c_ref_compare = 
      (&_category < &(e._category))
       ? (-1)
       : ((&_category > &(e._category)) ? 1 : 0);
   int c_compare = 
      (c_id.get() == NULL)
       ? ((e_c_id.get() == NULL) ? c_ref_compare : -1)
       : ((e_c_id.get() == NULL) ? 1 : c_id->compare_to(*e_c_id));
   /* compare exemplars */
   int e_compare = 
      (_id < e._id)
       ? (-1)
       : ((_id > e._id) ? 1 : 0);
   return (c_compare == 0) ? e_compare : c_compare;
}
   
} /* namespace ids */
} /* namespace models */
} /* namespace recognition */
} /* namespace vision */
