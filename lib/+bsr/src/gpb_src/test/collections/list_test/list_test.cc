/*
 * List test.
 */
#include "collections/abstract/collection.hh"
#include "collections/list.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"

#include <iostream>

using collections::abstract::collection;
using collections::list;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

int main() {
   /* basic list tests */
   list<const int> lst;
   lst.append(3).prepend(2).append(4).prepend(1);
   std::cout << lst.head() << " ";
   lst.remove_head();
   std::cout << lst.head() << " ";
   lst.remove_head();
   std::cout << lst.head() << " ";
   lst.remove_head();
   std::cout << lst.tail() << "\n";
   lst.remove_head();
   /* test list as collection */
   lst.append(3).prepend(2).append(4).prepend(1);
   collection<const int>& c = lst;
   c.add(c);
   /* test iterators over list */
   auto_ptr< iterator<const int> > i = c.iter_create();
   while (i->has_next()) {
      std::cout << i->next() << " ";
   }
   std::cout << "\n";
   for (list<const int>::iterator_reverse_t ii(lst); ii.has_next(); )
      std::cout << ii.next() << " ";
   std::cout << "\n";
   i = lst.iter_reverse_create();
   while (i->has_next()) {
      std::cout << i->next() << " ";
   }
   std::cout << "\n";
   /* test sorting list */
   lst.sort();
   for (list<const int>::iterator_t ii(lst); ii.has_next(); )
      std::cout << ii.next() << " ";
   std::cout << "\n";
   /* test sorting more complicated list */
   list<const double> a;
   a.append(0);
   a.append(0.166035);
   a.append(0.0381288);
   a.append(0.00558394);
   a.append(0.191116);
   a.append(0.0099273);
   a.append(0.198789);
   a.append(0.23788);
   a.append(0.1603);
   a.append(0.264449);
   a.append(0.297406);
   a.append(0.0491625);
   a.append(0.155613);
   a.append(0.346112);
   a.append(0.119747);
   a.append(0.340048);
   a.append(0.13701);
   a.append(0.314217);
   a.append(0);
   a.append(0);
   list<const double> b(a);
   a.sort();
   for (list<const double>::iterator_t ii(a); ii.has_next(); )
      std::cout << ii.next() << " ";
   std::cout << "\n";
   b.unique();
   for (list<const double>::iterator_t ii(b); ii.has_next(); )
      std::cout << ii.next() << " ";
   std::cout << "\n";
   return 0;
}
