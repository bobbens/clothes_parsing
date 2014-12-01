/*
 * Multimap test.
 */
#include "collections/abstract/collection.hh"
#include "collections/list.hh"
#include "collections/multimap.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"

#include <iostream>

using collections::abstract::collection;
using collections::list;
using collections::multimap;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

int main() {
   multimap<const int, const double> m;
   m.add(3, 1.23).add(2, 8.0).add(4, 156.3).add(1, 0.2).add(2, 6.7);
   auto_ptr< iterator<const int> > i = m.iter_create();
   while (i->has_next()) {
      const int& t = i->next();
      list<const double> im = m.find_images(t);
      auto_ptr< iterator<const double> > ii = im.iter_create();
      while (ii->has_next()) {
         const double& t_im = ii->next();
         std::cout << t << " -> " << t_im << "\n";
      }
   }
   auto_ptr< list<const int> >    dmn;
   auto_ptr< list<const double> > rng;
   m.contents(dmn, rng);
   auto_ptr< iterator<const int> > i_dmn = dmn->iter_create();
   auto_ptr< iterator<const double> > i_rng = rng->iter_create();
   while (i_dmn->has_next())
      std::cout << i_dmn->next() << " -> " << i_rng->next() << "\n";
   return 0;
}
