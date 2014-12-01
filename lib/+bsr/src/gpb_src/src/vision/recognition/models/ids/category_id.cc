/*
 * Category identity.
 */
#include "lang/pointers/safe_ptr.hh"
#include "vision/recognition/models/category.hh"
#include "vision/recognition/models/ids/category_id.hh"

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
category_id::category_id(category& c, unsigned long id)
 : _category(c),
   _id(id)
{ }

/*
 * Copy constructor.
 */
category_id::category_id(const category_id& c)
 : _category(c._category),
   _id(c._id)
{ }

/*
 * Destructor.
 */
category_id::~category_id() {
   /* do nothing */
}

/*
 * Get the parent category to which the subcategory belongs.
 */
category& category_id::parent_category() const {
   return _category;
}

/*
 * Get the subcategory number within the parent category.
 */
unsigned long category_id::id() const {
   return _id;
}

/*
 * Compare to another category identity.
 */   
int category_id::compare_to(const category_id& c) const {
   /* compare parent categories */
   safe_ptr<const category_id>   c_id =   _category.category_identity();
   safe_ptr<const category_id> c_c_id = c._category.category_identity();
   int parent_c_ref_compare = 
      (&_category < &(c._category))
       ? (-1)
       : ((&_category > &(c._category)) ? 1 : 0);
   int parent_c_compare = 
      (c_id.get() == NULL)
       ? ((c_c_id.get() == NULL) ? parent_c_ref_compare : -1)
       : ((c_c_id.get() == NULL) ? 1 : c_id->compare_to(*c_c_id));
   /* compare categories */
   int c_compare = 
      (_id < c._id)
       ? (-1)
       : ((_id > c._id) ? 1 : 0);
   return (parent_c_compare == 0) ? c_compare : parent_c_compare;
}

} /* namespace ids */
} /* namespace models */
} /* namespace recognition */
} /* namespace vision */
