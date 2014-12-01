/*
 * Hash map test.
 */
#include "collections/abstract/collection.hh"
#include "collections/hash_map.hh"
#include "functors/equalable_functors.hh"
#include "functors/comparable_functors.hh"
#include "functors/mappable_functors.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"

#include <iostream>

using collections::abstract::collection;
using collections::hash_map;
using functors::equal_functors;
using functors::compare_functors;
using functors::mappable_functor;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

class h_fun : public mappable_functor<const int, const int> {
public:
   using mappable_functor<const int, const int>::operator();
   const int& operator()(const int& i) const {
      return i;
   }
};

int main() {
   /* basic hash map tests */
   h_fun h_func;
   hash_map<const int, const double, const int> m(
      equal_functors<const int>::f_equal(),
      compare_functors<const int>::f_compare(),
      h_func);
   m.add(3, 1.23).add(2, 8.0).add(4, 6.7).add(1, 0.2);
   auto_ptr< iterator<const int> > i = m.iter_create();
   while (i->has_next()) {
      const int& t = i->next();
      std::cout << t << " -> " << m.find_image(t) << "\n";
   }
   return 0;
}
