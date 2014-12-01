/*
 * Hash set test.
 */
#include "collections/abstract/collection.hh"
#include "collections/hash_set.hh"
#include "functors/equalable_functors.hh"
#include "functors/comparable_functors.hh"
#include "functors/mappable_functors.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"

#include <iostream>

using collections::abstract::collection;
using collections::hash_set;
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
   /* basic hash set tests */
   h_fun h_func;
   hash_set<const int, const int> s(
      equal_functors<const int>::f_equal(),
      compare_functors<const int>::f_compare(),
      h_func);
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
   for (hash_set<const int, const int>::iterator_reverse_t ii(s); ii.has_next(); )
      std::cout << ii.next() << " ";
   std::cout << "\n";
   i = s.iter_reverse_create();
   while (i->has_next()) {
      std::cout << i->next() << " ";
   }
   std::cout << "\n";
   return 0;
}
