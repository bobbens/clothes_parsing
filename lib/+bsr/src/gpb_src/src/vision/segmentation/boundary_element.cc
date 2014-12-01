/*
 * Boundary element.
 */
#include "collections/array_list.hh"
#include "collections/pointers/auto_collection.hh"
#include "lang/array.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/geometry/point_2D.hh"
#include "math/libraries/lib_image.hh"
#include "math/math.hh"
#include "math/matrices/matrix.hh"
#include "vision/segmentation/boundary_element.hh"

namespace vision {
namespace segmentation {
/*
 * Imports.
 */
using collections::array_list;
using collections::pointers::auto_collection;
using lang::array;
using lang::pointers::auto_ptr;
using math::geometry::point_2D;
using math::libraries::lib_image;
using math::matrices::matrix;

/*
 * Create boundary elements for the given collection of contours.
 */
auto_collection< boundary_element, array_list<boundary_element> > 
   boundary_element::create_boundary_elements(
      const matrix<unsigned long>&                 contour_assignments,
      const array_list<lib_image::contour_vertex>& vertices,
      const array_list<lib_image::contour_edge>&   edges)
{
   /* get number of vertices/edges */
   unsigned long n_vertices = vertices.size();
   unsigned long n_edges = edges.size();
   /* count the number of constraint edges out of each vertex */
   array<unsigned long> v_constraint_count(n_vertices);
   for (unsigned long n = 0; n < n_edges; n++) {
      const lib_image::contour_edge& e = edges[n];
      if (!(e.is_completion)) {
         v_constraint_count[e.vertex_start->id]++;
         v_constraint_count[e.vertex_end->id]++;
      }
   }
   /* initialize collection of boundary elements */
   auto_collection< boundary_element, array_list<boundary_element> >
      boundary_elements(new array_list<boundary_element>());
   /* create boundary element for each contour */
   unsigned long size_x = contour_assignments.size(0);
   unsigned long size_y = contour_assignments.size(1);
   for (unsigned long n = 0; n < n_edges; n++) {
      /* create boundary element */
      const lib_image::contour_edge& e = edges[n];
      auto_ptr<boundary_element> b(new boundary_element(e));
      /* set flags */
      unsigned long n_constraint_start = v_constraint_count[e.vertex_start->id];
      unsigned long n_constraint_end   = v_constraint_count[e.vertex_end->id];
      unsigned long e_constraint = e.is_completion ? 0 : 1;
      b->_is_intruder =
         (e.is_completion) && 
         ((n_constraint_start == 2) || (n_constraint_end == 2));
      b->_is_interior =
         ((n_constraint_start - e_constraint) > 0) &&
         ((n_constraint_end - e_constraint) > 0);
      b->_is_border = true;
      unsigned long n_pts = e.x_coords.size();
      for (unsigned long p = 0; (p < n_pts) && (b->_is_border); p++) {
         bool border_x = (e.x_coords[p] == 0) || ((e.x_coords[p]+1) == size_x);
         bool border_y = (e.y_coords[p] == 0) || ((e.y_coords[p]+1) == size_y);
         if (!(border_x || border_y))
            b->_is_border = false;
      }
      /* add boundary element to set */
      boundary_elements->add(*b);
      b.release();
   }
   return boundary_elements;
}

/*
 * Protected constructor.
 * Create a boundary element corresponding to the given contour edge.
 * Initialize its status flags to false.
 */
boundary_element::boundary_element(const lib_image::contour_edge& e)
 : _edge(e),
   _is_intruder(false),
   _is_interior(false),
   _is_border(false)
{ }
 
/*
 * Copy constructor.
 */
boundary_element::boundary_element(const boundary_element& b)
 : _edge(b._edge),
   _is_intruder(b._is_intruder),
   _is_interior(b._is_interior),
   _is_border(b._is_border)
{ }

/*
 * Destructor.
 */
boundary_element::~boundary_element() {
   /* do nothing */
}

/*
 * Get the id of the smooth contour represented by the boundary element.
 */
unsigned long boundary_element::contour_id() const {
   return _edge.id;
}

/*
 * Get the id of the contour set to which the boundary element belongs.
 *
 * Each contour set is a longest possible sequence of connected contours 
 * such that each connection vertex joins two constraint boundaries and 
 * all other boundaries intersecting it are completion edges.
 *
 * Contour sets correspond to the image contours that existed prior to
 * subdivision into approximate straight line segments.  Each completion
 * edge is in its own separate contour set.
 */
unsigned long boundary_element::contour_set_id() const {
   return _edge.contour_equiv_id;
}

/*
 * Check whether the boundary element is a completion edge.
 */
bool boundary_element::is_completion() const {
   return _edge.is_completion;
}

/*
 * Check whether the boundary element is an intruding edge.
 * An intruding edge is a completion edge which interrupts a contour set.
 */
bool boundary_element::is_intruder() const {
   return _is_intruder;
}

/*
 * Check whether the boundary element is an interior edge to two constraint
 * edges.  An interior edge is connected to at least one constraint edge at
 * each of its vertices.
 */
bool boundary_element::is_interior() const {
   return _is_interior;
}

/*
 * Check whether the boundary element lies entirely along the image border.
 */
bool boundary_element::is_border() const {
   return _is_border;
}

/*
 * Get the number of pixels along the boundary element's smooth contour.
 */
unsigned long boundary_element::size() const {
   return _edge.x_coords.size();
}

/*
 * Get the length of the straight line segment approxmation to the boundary.
 */
double boundary_element::length() const {
   unsigned long n_pts = _edge.x_coords.size();
   point_2D p_start(
      static_cast<double>(_edge.x_coords[0]),
      static_cast<double>(_edge.y_coords[0])
   );
   point_2D p_end(
      static_cast<double>(_edge.x_coords[n_pts-1]),
      static_cast<double>(_edge.y_coords[n_pts-1])
   );
   return abs(p_end - p_start);
}

/*
 * Get the orientation of the normal to the boundary.
 * The returned orientation is in the range [0, pi).
 */
double boundary_element::orientation() const {
   unsigned long n_pts = _edge.x_coords.size();
   point_2D p_start(
      static_cast<double>(_edge.x_coords[0]),
      static_cast<double>(_edge.y_coords[0])
   );
   point_2D p_end(
      static_cast<double>(_edge.x_coords[n_pts-1]),
      static_cast<double>(_edge.y_coords[n_pts-1])
   );
   point_2D p_diff = p_end - p_start;
   double ang = math::atan2(p_diff.y(), p_diff.x());
   if (ang < 0)
      ang += M_PI;
   return ang;
}

/*
 * Get x-coordinates of points along the boundary element.
 */
const array<unsigned long>& boundary_element::x_coords() const {
   return _edge.x_coords;
}

/*
 * Get y-coordinates of points along the boundary element.
 */
const array<unsigned long>& boundary_element::y_coords() const {
   return _edge.y_coords;
}

} /* namespace segmentation */
} /* namespace vision */
