/*
 * Line segment intersection test.
 */
#include "collections/list.hh"
#include "io/streams/cout.hh"
#include "lang/array.hh"
#include "lang/exceptions/exception.hh"
#include "lang/exceptions/ex_index_out_of_bounds.hh"
#include "math/geometry/point_2D.hh"
#include "math/geometry/seg_intersect.hh"

using collections::list;
using io::streams::cout;
using lang::array;
using lang::exceptions::exception;
using lang::exceptions::ex_index_out_of_bounds;
using math::geometry::point_2D;
using math::geometry::seg_intersect;

int main() {
   try {
      /* create point set */
      point_2D p0(0,0);
      point_2D p1(0,1);
      point_2D p2(1,0);
      point_2D p3(1,1);
      point_2D p4(1.5,0.5);
      point_2D p5(0,-1);
      point_2D p6(0,2);
      point_2D p7(0,3);
      point_2D p8(3,1);
      list<point_2D> points;
      points.add(p0); 
      points.add(p1);
      points.add(p2);
      points.add(p3);
      points.add(p4);
      points.add(p5);
      points.add(p6);
      points.add(p7);
      points.add(p8);
      /* create segments */
      array<unsigned long> seg_start(12);
      array<unsigned long> seg_end(12);
      seg_start[0] = 0; seg_end[0] = 1;
      seg_start[1] = 0; seg_end[1] = 2;
      seg_start[2] = 0; seg_end[2] = 3;
      seg_start[3] = 7; seg_end[3] = 3;
      seg_start[4] = 1; seg_end[4] = 3;
      seg_start[5] = 2; seg_end[5] = 3;
      seg_start[6] = 0; seg_end[6] = 4;
      seg_start[7] = 5; seg_end[7] = 6;
      seg_start[8] = 5; seg_end[8] = 8;
      seg_start[9] = 1; seg_end[9] = 2;
      seg_start[10] = 5; seg_end[10] = 7;
      seg_start[11] = 0; seg_end[11] = 8;
      /* compute intersections */
      seg_intersect s(points, seg_start, seg_end);
      /* output intersections */
      unsigned long n_v = s.vertices_size();
      for (unsigned long n = 0; n < n_v; n++) {
         cout << s.vertex(n) << "\n";
         cout << "   endpoint: ";
         array<unsigned long> s_ids = s.endpoint_intersection(n);
         for (unsigned long j = 0; j < s_ids.size(); j++) {
            unsigned long id = s_ids[j];
            cout << "{" 
                 << s.vertex(s.segment_vertex_id(id,0)) << ", "
                 << s.vertex(s.segment_vertex_id(id,1)) << "}"
                 << "  ";
         }
         cout << "\n";
         cout << "   interior: ";
         s_ids = s.interior_intersection(n);
         for (unsigned long j = 0; j < s_ids.size(); j++) {
            unsigned long id = s_ids[j];
            cout << "{" 
                 << s.vertex(s.segment_vertex_id(id,0)) << ", "
                 << s.vertex(s.segment_vertex_id(id,1)) << "}"
                 << "  ";
         }
         cout << "\n";
      }
   } catch (ex_index_out_of_bounds& e) {
      cout << e.index() << "\n";
      cout << e << "\n";
   } catch (exception& e) {
      cout << e << "\n";
   }
   return 0;
}
