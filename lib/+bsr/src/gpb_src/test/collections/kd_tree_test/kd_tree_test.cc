/*
 * Kd tree test.
 */
#include "collections/abstract/collection.hh"
#include "collections/list.hh"
#include "collections/queue.hh"
#include "collections/kd_tree.hh"
#include "interfaces/comparable.hh"
#include "interfaces/keyable.hh"
#include "interfaces/kd_treeable.hh"
#include "lang/exceptions/exception.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/math.hh"

#include <iostream>

using collections::abstract::collection;
using collections::list;
using collections::queue;
using collections::kd_tree;
using interfaces::comparable;
using interfaces::keyable;
using interfaces::kd_tree_key;
using interfaces::kd_treeable;
using lang::exceptions::exception;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

/*
 * A kd-treeable point class.
 */
class point : public comparable<point>,
              public kd_treeable<point> {
public:
   /*
    * Constructor.
    */
   point(double x, double y);

   /*
    * Copy constructor.
    */
   point(const point&);

   /*
    * Destructor.
    */
   virtual ~point();

   /*
    * Distance from point to key.
    */
   double distance_to(const kd_tree_key<>&) const;

   /*
    * Distance from point to another point.
    */
   double distance_to(const point&) const;

   /*
    * Compare to a key.
    */
   int compare_to(const kd_tree_key<>&) const;

   /*
    * Lexicographically compare to another point.
    */
   int compare_to(const point&) const;

   /*
    * Print point coordinates.
    */ 
   void print() const;
  
   /*
    * Dimensionality of points.
    */
   static unsigned long dimensionality(const collection<point>&);

   /*
    * Compute key for split collection of points.
    */
   static kd_tree_key<> key_split(const collection<point>&);

protected:
   double x;
   double y;
};

/*
 * Constructor.
 */
point::point(double x, double y) : x(x), y(y) { }

/*
 * Copy constructor.
 */
point::point(const point& p) : x(p.x), y(p.y) { }

/*
 * Destructor.
 */
point::~point() {
   /* do nothing */
}
 
/*
 * Distance from point to key.
 */
double point::distance_to(const kd_tree_key<>& k) const {
   double dist = 0;
   if (k.split_dimension() == 0)
      dist = this->x - k.split_value();
   else
      dist = this->y - k.split_value();
   return ((dist < 0) ? (-dist) : dist);
}
   
/*
 * Distance from point to another point.
 */
double point::distance_to(const point& p) const {
   double dist_x = this->x - p.x;
   double dist_y = this->y - p.y;
   return math::sqrt(dist_x*dist_x + dist_y*dist_y);
}

/*
 * Compare to a key.
 */
int point::compare_to(const kd_tree_key<>& k) const {
   double dist = 0;
   if (k.split_dimension() == 0)
      dist = this->x - k.split_value();
   else
      dist = this->y - k.split_value();
   return ((dist < 0) ? (-1) : 
           ((dist > 0) ? 1 : 0));
}

/*
 * Lexicographically compare to another point.
 */
int point::compare_to(const point& p) const {
   if (this->x < p.x) {
      return (-1);
   } else if (this->x > p.x) {
      return 1;
   } else {
      return ((this->y < p.y) ? (-1) : 
              ((this->y > p.y) ? 1 : 0));
   }
}

/*
 * Print point coordinates.
 */
void point::print() const {
   std::cout << "(" << this->x << ", " << this->y << ")";
}

/*
 * Points are two-dimensional.
 */
unsigned long point::dimensionality(const collection<point>&) {
   return 2;
}

/*
 * Compute key for split collection of points.
 */
kd_tree_key<> point::key_split(const collection<point>& c) {
   unsigned long split_dim;
   double split_val;
   queue<point> q(c);
   unsigned long size = q.size();
   for (unsigned long i = 0; i < (size/2)-1; i++)
      q.dequeue();
   point& p_left = q.dequeue();
   point& p_right = q.dequeue();
   if (p_left.x != p_right.x) {
      split_dim = 0;
      split_val = p_right.x;
   } else {
      split_dim = 1;
      split_val = p_right.y;
   }
   return kd_tree_key<>(split_dim, split_val);
}

/*
 * Print a collection of points.
 */
void print_collection(const collection<point>& c) {
   auto_ptr< iterator<point> > i = c.iter_create();
   while (i->has_next()) {
      i->next().print();
      std::cout << " ";
   }
   std::cout << "\n";
}

/*
 * Delete a collection of points.
 */
void delete_collection(const collection<point>& c) {
   auto_ptr< iterator<point> > i = c.iter_create();
   while (i->has_next()) {
      point& p = i->next();
      delete &p;
   }
}

/*
 * Test kd-tree of points.
 */
int main() {
   /* build a kd-tree */
   list<point> lst;
   lst
      .append(*(new point(1,1)))
      .append(*(new point(2,3)))
      .append(*(new point(-4,3.2)))
      .append(*(new point(10,12)))
      .append(*(new point(10,11)))
      .append(*(new point(-0.5,1)));
   std::cout << "list = ";
   print_collection(lst);
   kd_tree<point> k_tree(lst);
   std::cout << "kd tree = ";
   print_collection(k_tree);
   /* test finding a point */
   point p1(-4,3.2);
   point& p1_found = k_tree.find(p1);
   std::cout << "found ";
   p1_found.print();
   std::cout << "\n";
   /* test nearest neighbors */
   point p2(-2.1,3.49);
   point& p2_nn = k_tree.find_nn(p2);
   p2.print();
   std::cout << " nearest's neighbor is ";
   p2_nn.print();
   std::cout << "\n";
   delete_collection(lst);
   return 0;
}
