/*
 * Array list test.
 */
#include "collections/abstract/collection.hh"
#include "collections/array_list.hh"
#include "lang/array.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"

#include <iostream>

using collections::abstract::collection;
using collections::array_list;
using lang::array;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

int main() {
   /* basic array list tests */
   array_list<const int> lst;
   lst.append(3).prepend(2).append(4).prepend(1);
   std::cout << lst.head() << " ";
   lst.remove_head();
   std::cout << lst.head() << " ";
   lst.remove_head();
   std::cout << lst.head() << " ";
   lst.remove_head();
   std::cout << lst.tail() << "\n";
   lst.remove_head();
   /* test array list as collection */
   lst.append(3).prepend(2).append(4).prepend(1);
   collection<const int>& c = lst;
   c.add(c);
   /* test iterators over array list */
   auto_ptr< iterator<const int> > i = c.iter_create();
   while (i->has_next()) {
      std::cout << i->next() << " ";
   }
   std::cout << "\n";
   for (array_list<const int>::iterator_reverse_t ii(lst); ii.has_next(); )
      std::cout << ii.next() << " ";
   std::cout << "\n";
   i = lst.iter_reverse_create();
   while (i->has_next()) {
      std::cout << i->next() << " ";
   }
   std::cout << "\n";
   /* test sorting array list */
   {
      collections::abstract::array<const int>& a(lst);
      a.sort();
      for (array_list<const int>::iterator_t ii(lst); ii.has_next(); )
         std::cout << ii.next() << " ";
      std::cout << "\n";
   }
   /* test uniqueness */
   {
      array_list<const int> a(c);
      a.reverse();
      a.unique();
      for (array_list<const int>::iterator_t ii(a); ii.has_next(); )
         std::cout << ii.next() << " ";
      std::cout << "\n";
   }
   /* test unique_idx */
   {
      array_list<const int> a(c);
      a.reverse();
      array<unsigned long> idx = a.unique_idx();
      for (array_list<const int>::iterator_t ii(a); ii.has_next(); )
         std::cout << ii.next() << " ";
      std::cout << "\n";
      for (unsigned long n = 0; n < idx.size(); n++)
         std::cout << idx[n] << " ";
      std::cout << "\n";
   }
   return 0;
}
