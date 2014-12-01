/*
 * Queue test.
 */
#include "collections/queue.hh"
#include "functors/comparable_functors.hh"
#include "lang/exceptions/exception.hh"
#include <iostream>

using collections::queue;
using functors::compare_functors;
using lang::exceptions::exception;

queue<const int>& create_queue() {
   return *(new queue<const int>(compare_functors<const int>::f_compare_reverse()));
}

int main() {
   try {
      queue<const int>& q = create_queue();
      q.enqueue(3).enqueue(2).enqueue(4).enqueue(1);
      std::cout << q.dequeue() << "\n";
      std::cout << q.dequeue() << "\n";
      std::cout << q.dequeue() << "\n";
      std::cout << q.dequeue() << "\n";
      delete &q;
   } catch (exception& e) {
      std::cout << e << "\n";
      throw;
   }
   return 0;
}
