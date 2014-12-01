/*
 * Set test.
 */
#include "collections/abstract/collection.hh"
#include "collections/set.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"

#include <iostream>

using collections::abstract::collection;
using collections::set;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

int main() {
   /* basic splay set tests */
   set<const int> s;
   s.add(3).add(2).add(4).add(1);
   /* test set as collection */
   collection<const int>& c = s;
   c.add(c);
   /* test iterators over set */
   auto_ptr< iterator<const int> > i = c.iter_create();
   while (i->has_next()) {
      std::cout << i->next() << " ";
   }
   std::cout << "\n";
   for (set<const int>::iterator_reverse_t ii(s); ii.has_next(); )
      std::cout << ii.next() << " ";
   std::cout << "\n";
   i = s.iter_reverse_create();
   while (i->has_next()) {
      std::cout << i->next() << " ";
   }
   std::cout << "\n";
   return 0;
}
