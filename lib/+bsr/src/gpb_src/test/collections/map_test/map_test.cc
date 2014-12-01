/*
 * Map test.
 */
#include "collections/abstract/collection.hh"
#include "collections/map.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"

#include <iostream>

using collections::abstract::collection;
using collections::map;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

int main() {
   map<const int, const double> m;
   m.add(3, 1.23).add(2, 8.0).add(4, 6.7).add(1, 0.2);
   auto_ptr< iterator<const int> > i = m.iter_create();
   while (i->has_next()) {
      const int& t = i->next();
      std::cout << t << " -> " << m.find_image(t) << "\n";
   }
   return 0;
}
