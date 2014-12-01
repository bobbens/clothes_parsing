/*
 * Triangulation test.
 */
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "io/streams/cin.hh"
#include "io/streams/cout.hh"
#include "lang/array.hh"
#include "lang/exceptions/exception.hh"
#include "lang/exceptions/ex_index_out_of_bounds.hh"
#include "lang/null.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/geometry/point_2D.hh"
#include "math/geometry/triangulation.hh"

#include <cstdlib>

using collections::list;
using collections::pointers::auto_collection;
using io::streams::cin;
using io::streams::cout;
using lang::array;
using lang::exceptions::exception;
using lang::exceptions::ex_index_out_of_bounds;
using lang::pointers::auto_ptr;
using math::geometry::point_2D;
using math::geometry::triangulation;

int main() {
   try {
      /* load point set */
      auto_collection< point_2D, list<point_2D> > points(new list<point_2D>());
      char line[100];
      cin.getline(line, 100);
      unsigned long n_points = strtol(line,NULL,10);
      for (unsigned long n = 0; n < n_points; n++) {
         cin.getline(line,100);
         char* curr = line;
         strtol(curr,&curr,10);
         double x = strtod(curr,&curr);
         double y = strtod(curr,&curr);
         auto_ptr<point_2D> p(new point_2D(x,y));
         points->add(*p);
         p.release();
      }
      auto_ptr<triangulation> t_map = triangulation::delaunay(*points);
      unsigned long n_t = t_map->triangles_size();
      cout << n_t << " 3 0\n";
      for (unsigned long n = 0; n < n_t; n++) {
         cout << n+1 << " "
              << t_map->triangle_vertex_id(n,0)+1 << " " 
              << t_map->triangle_vertex_id(n,1)+1 << " " 
              << t_map->triangle_vertex_id(n,2)+1 << "\n" ;
      }
   } catch (ex_index_out_of_bounds& e) {
      cout << e.index() << "\n";
      cout << e << "\n";
   } catch (exception& e) {
      cout << e << "\n";
   }
   return 0;
}
