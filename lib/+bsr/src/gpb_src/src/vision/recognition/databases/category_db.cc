/*
 * Category DB.
 */
#include "collections/list.hh"
#include "collections/kd_tree.hh"
#include "lang/iterators/iterator.hh"
#include "lang/null.hh"
#include "lang/pointers/auto_ptr.hh"
#include "vision/features/geometric_blur.hh"
#include "vision/recognition/databases/category_db.hh"
#include "vision/recognition/models/category.hh"
#include "vision/recognition/models/exemplar.hh"

namespace vision {
namespace recognition {
namespace databases {
/*
 * Imports.
 */
using collections::list;
using collections::kd_tree;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;
using vision::features::geometric_blur;
using vision::recognition::models::category;
using vision::recognition::models::exemplar;


/*
 * Constructor.
 * Build a database for the given category.
 */
category_db::category_db(auto_ptr<category> c)
 : _root_category(c),
   _gb_db(NULL)
{
   /* build database search structures */
   this->build_db();
}

/*
 * Copy constructor.
 * Make a deep copy of the database.
 */
category_db::category_db(const category_db& c_db)
 : _root_category(new category(*(c_db._root_category))),
   _gb_db(NULL)
{
   /* build database search structures */
   this->build_db();
}

/*
 * Build the database search structures from the root category.
 */
void category_db::build_db() {
   /* initialize list of categories to process */
   list<category> categories;
   categories.add(*_root_category);
   /* initialize lists of all features in database */
   list<geometric_blur> gb_features;
   /* collect features from each exemplar */
   while (!(categories.is_empty())) {
      /* grab next category */
      category& c = categories.remove_head();
      categories.add(c.subcategories());
      /* process exemplars in category */
      auto_ptr< iterator<exemplar> > i = c.exemplars().iter_create();
      while (i->has_next()) {
         exemplar& e = i->next();
         gb_features.add(e.gb_features());
      }
   }
   /* build search structures */
   _gb_db.reset(new kd_tree<geometric_blur>(gb_features));
}

/*
 * Destructor.
 */
category_db::~category_db() {
   /* do nothing */
}

/*
 * Database search structure access.
 */
const kd_tree<geometric_blur>& category_db::gb_features() const {
   return *_gb_db;
}

} /* namespace databases */
} /* namespace recognition */
} /* namespace vision */
