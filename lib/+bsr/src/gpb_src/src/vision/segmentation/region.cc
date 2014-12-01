/*
 * Region.
 */
#include "vision/segmentation/region.hh"

namespace vision {
namespace segmentation {

/*
 * Constructor.
 */
region::region()
 : _id(0),
   _boundary_map(),
   _size(0),
   _sum_L(0),
   _sum_L2(0),
   _sum_a(0),
   _sum_a2(0),
   _sum_b(0),
   _sum_b2(0)
{ }

/*
 * Copy constructor.
 */
region::region(const region& r)
 : _id(r._id),
   _boundary_map(r._boundary_map),
   _size(r._size),
   _sum_L(r._sum_L),
   _sum_L2(r._sum_L2),
   _sum_a(r._sum_a),
   _sum_a2(r._sum_a2),
   _sum_b(r._sum_b),
   _sum_b2(r._sum_b2)
{ }

/*
 * Destructor.
 */
region::~region() {
   /* do nothing */
}
 
} /* namespace segmentation */
} /* namespace vision */
