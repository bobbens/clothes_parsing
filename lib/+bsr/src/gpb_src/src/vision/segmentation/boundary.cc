/*
 * Boundary.
 */
#include "vision/segmentation/boundary.hh"

namespace vision {
namespace segmentation {

/*
 * Constructor.
 */
boundary::boundary()
 : _id(0),
   _size(0),
   _sum_contrast(0),
   _sum_pb(0)
{ }

/*
 * Destructor.
 */
boundary::~boundary() {
   /* do nothing */
}

} /* namespace segmentation */
} /* namespace vision */
