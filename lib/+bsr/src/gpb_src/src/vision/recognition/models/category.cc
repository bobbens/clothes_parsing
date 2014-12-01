/*
 * Category.
 */
#include "collections/array_list.hh"
#include "collections/pointers/auto_collection.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"
#include "lang/pointers/safe_ptr.hh"
#include "vision/recognition/models/exemplar.hh"
#include "vision/recognition/models/category.hh"
#include "vision/recognition/models/ids/category_id.hh"

namespace vision {
namespace recognition {
namespace models {
/*
 * Imports.
 */
using collections::array_list;
using collections::pointers::auto_collection;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;
using lang::pointers::safe_ptr;
using vision::recognition::models::ids::category_id;

/*
 * Constructor.
 * Create an empty category and set its identity.
 */
category::category(auto_ptr<category_id> c_id)
 : _exemplars(new array_list<exemplar>),
   _subcategories(new array_list<category>),
   _c_id(c_id)
{ }

/*
 * Copy constructor.
 * Perform a deep copy of the category.
 * Assign the copy the specified id.
 */
category::category(const category& c, auto_ptr<category_id> c_id)
 : _exemplars(new array_list<exemplar>),
   _subcategories(new array_list<category>),
   _c_id(c_id)
{
   /* copy exemplars */
   unsigned long n = 0;
   for (auto_ptr< iterator<exemplar> > i = c._exemplars->iter_create(); i->has_next(); n++) {
      auto_ptr<exemplar_id> e_id(new exemplar_id(*this, n));
      auto_ptr<exemplar> e_copy(new exemplar(i->next(), e_id));
      _exemplars->add(*e_copy);
      e_copy.release();
   }  
   /* copy subcategories */
   n = 0;
   for (auto_ptr< iterator<category> > i = c._subcategories->iter_create(); i->has_next(); n++) {
      auto_ptr<category_id> c_id(new category_id(*this, n));
      auto_ptr<category> c_copy(new category(i->next(), c_id));
      _subcategories->add(*c_copy);
      c_copy.release();
   }
}

/*
 * Destructor.
 */
category::~category() {
   /* do nothing */
}

/*
 * Add an exemplar to the category.
 * Return a reference to the category.
 */
category& category::add_exemplar(auto_ptr<exemplar> e) {
   auto_ptr<exemplar_id> e_id(new exemplar_id(*this, _exemplars->size()));
   e->exemplar_identity(e_id);
   _exemplars->add(*e);
   e.release();
   return *this;
}

/*
 * Add a subcategory to the category.
 * Return a reference to the category.
 */
category& category::add_subcategory(auto_ptr<category> c) {
   auto_ptr<category_id> c_id(new category_id(*this, _subcategories->size()));
   c->category_identity(c_id);
   _subcategories->add(*c);
   c.release();
   return *this;
}

/*
 * Get exemplars within the category.
 */
const array_list<exemplar>& category::exemplars() const {
   return *_exemplars;
}

/*
 * Get subcategories.
 */
const array_list<category>& category::subcategories() const {
   return *_subcategories;
}
 
/*
 * Get category identity.
 */
safe_ptr<const category_id> category::category_identity() const {
   return safe_ptr<const category_id>(_c_id.get());
}

/*
 * Set category identity.
 * Return the new identity.
 */
safe_ptr<const category_id> category::category_identity(auto_ptr<category_id> c_id) {
   _c_id = c_id;
   return safe_ptr<const category_id>(_c_id.get());
}
   
} /* namespace models */
} /* namespace recognition */
} /* namespace vision */
