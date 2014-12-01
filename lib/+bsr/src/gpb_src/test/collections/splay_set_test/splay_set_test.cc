/*
 * Splay set test.
 */
#include "collections/abstract/collection.hh"
#include "collections/splay_set.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"

#include <iostream>

using collections::abstract::collection;
using collections::splay_set;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

int main() {
   /* basic splay set tests */
   splay_set<const int> s;
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
   for (splay_set<const int>::iterator_reverse_t ii(s); ii.has_next(); )
      std::cout << ii.next() << " ";
   std::cout << "\n";
   i = s.iter_reverse_create();
   while (i->has_next()) {
      std::cout << i->next() << " ";
   }
   std::cout << "\n";
   /* test finding min/max element */
   std::cout << s.find_min() << "\n";
   std::cout << s.find_max() << "\n";
   /* test containst prev/next */
   std::cout << (s.contains_prev(3) ? "true" : "false") << "\n";
   std::cout << (s.contains_prev(1) ? "true" : "false") << "\n";
   std::cout << (s.contains_prev(4) ? "true" : "false") << "\n";
   std::cout << (s.contains_next(4) ? "true" : "false") << "\n";
   /* test finding prev/next element */
   std::cout << s.find_prev(3) << "\n";
   std::cout << s.find_next(2) << "\n";
   return 0;
}
