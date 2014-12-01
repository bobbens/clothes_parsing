/*
 * Key tree test.
 */
#include "collections/abstract/collection.hh"
#include "collections/list.hh"
#include "collections/queue.hh"
#include "collections/key_tree.hh"
#include "functors/equalable_functors.hh"
#include "functors/keyable_functors.hh"
#include "lang/exceptions/exception.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"

#include <iostream>

using collections::abstract::collection;
using collections::list;
using collections::queue;
using collections::key_tree;
using functors::equal_functor;
using functors::keyable_compare_functor;
using functors::keyable_split_functor;
using lang::exceptions::exception;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

/*
 * Specialize the keyable functors to integers.
 */
class my_compare_functor : public keyable_compare_functor<const int, int> {
public:
   int operator()(const int& i0, const int& i1) const {
      return i0-i1;
   }
};

class my_split_functor : public keyable_split_functor<const int, int> {
public:
   int operator()(const collection<const int>& c) const {
      queue<const int> q(c);
      unsigned long size = q.size();
      for (unsigned long i = 0; i < (size/2); i++)
         q.dequeue();
      const int& j = q.dequeue();
      while (!(q.is_empty()))
         q.dequeue();
      return j;
   }
};

void print_collection(const collection<const int>& c) {
   auto_ptr< iterator<const int> > i = c.iter_create();
   while (i->has_next()) 
      std::cout << i->next() << " ";
   std::cout << "\n";
}

int main() {
   list<const int> lst;
   lst.append(10).append(20).append(80).append(40).append(50).append(70).append(90).append(5);
   std::cout << "list = ";
   print_collection(lst);
   equal_functor<const int> eq_func;
   my_compare_functor my_c_func;
   my_split_functor my_s_func;
   key_tree<const int, int> k_tree(
      lst,
      eq_func,
      my_c_func,
      my_s_func
   );
   std::cout << "key tree = ";
   print_collection(k_tree);
   std::cout << "searching key tree... ";
   std::cout << k_tree.find(5) << " ";
   std::cout << k_tree.find(50) << " ";
   std::cout << k_tree.find(10) << " ";
   std::cout << k_tree.find(20) << " ";
   std::cout << k_tree.find(80) << " ";
   std::cout << k_tree.find(40) << " ";
   std::cout << k_tree.find(70) << " ";
   std::cout << k_tree.find(90) << " ";
   try {
      std::cout << k_tree.find(100) << "\n";
   } catch (exception& e) {
      std::cout << "\n" << e.what() << ": " << 100 << "\n";
   }
   return 0;
}
