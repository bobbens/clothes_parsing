/*
 * Feature.
 */
#include "lang/pointers/auto_ptr.hh"
#include "lang/pointers/safe_ptr.hh"
#include "vision/features/feature.hh"
#include "vision/features/feature_id.hh"

namespace vision {
namespace features {
/*
 * Imports.
 */
using lang::pointers::auto_ptr;
using lang::pointers::safe_ptr;

/*
 * Constructor.
 * Create a feature at the given location in scale space.
 * Optionally set its identity.
 */
feature::feature(
   double x_coord, 
   double y_coord, 
   double ori, 
   double scl,
   auto_ptr<feature_id> f_id)
 : _x(x_coord),
   _y(y_coord),
   _ori(ori),
   _scale(scl),
   _f_id(f_id)
{ }

/*
 * Copy constructor.
 * Assign the copy the specified id.
 */
feature::feature(
   const feature& f,
   auto_ptr<feature_id> f_id)
 : _x(f._x),
   _y(f._y),
   _ori(f._ori),
   _scale(f._scale),
   _f_id(f_id)
{ }

/*
 * Destructor.
 */
feature::~feature() {
   /* do nothing */
}

/*
 * Get x-coordinate of feature in image.
 */
double feature::x() const {
   return _x;
}

/*
 * Get y-coordinate of feature in image.
 */
double feature::y() const {
   return _y;
}

/*
 * Get local orientation (in radians) of feature.
 */
double feature::orientation() const {
   return _ori;
}

/*
 * Get scale of feature.
 */
double feature::scale() const {
   return _scale;
}

/*
 * Get feature identity.
 * Returns NULL if no identity has been assigned.
 */
safe_ptr<const feature_id> feature::feature_identity() const {
   return safe_ptr<const feature_id>(_f_id.get());
}

/*
 * Set feature identity.
 * Return the new identity.
 */
safe_ptr<const feature_id> feature::feature_identity(auto_ptr<feature_id> f_id) {
   _f_id = f_id;
   return safe_ptr<const feature_id>(_f_id.get());
}

} /* namespace features */
} /* namespace vision */
