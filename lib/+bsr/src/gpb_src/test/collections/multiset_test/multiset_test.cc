/*
 * Multiset test.
 */
#include "collections/abstract/collection.hh"
#include "collections/multiset.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"

#include <iostream>

using collections::abstract::collection;
using collections::multiset;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

int main() {
   /* basic multiset tests */
   multiset<const int> s;
   s.add(3).add(2).add(4).add(1).add(2).add(4);
   /* test multiset as collection */
   collection<const int>& c = s;
   c.add(c);
   /* test iterators over multiset */
   auto_ptr< iterator<const int> > i = c.iter_create();
   while (i->has_next()) {
      std::cout << i->next() << " ";
   }
   std::cout << "\n";
   for (multiset<const int>::iterator_reverse_t ii(s); ii.has_next(); )
      std::cout << ii.next() << " ";
   std::cout << "\n";
   i = s.iter_reverse_create();
   while (i->has_next()) {
      std::cout << i->next() << " ";
   }
   std::cout << "\n";
   return 0;
}
