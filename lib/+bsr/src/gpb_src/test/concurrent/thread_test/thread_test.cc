/*
 * Thread test.
 */
#include "collections/queue.hh"
#include "concurrent/threads/synchronization/synchronizables/synchronized.hh"
#include "concurrent/threads/runnable.hh"
#include "concurrent/threads/thread.hh"
#include "functors/comparable_functors.hh"
#include "interfaces/comparable.hh"

#include <iostream>

using collections::queue;
using concurrent::threads::synchronization::synchronizables::synchronized;
using concurrent::threads::runnable;
using concurrent::threads::thread;
using functors::compare_functors;
using interfaces::comparable;

class item : public comparable<item> {
public:
   item(int i) {
      this->i = i;
   }
   
   virtual int compare_to(const item& it) const {
      return (this->i - it.i);
   }
   
   void print() {
      std::cout << i << " ";
   }
protected:
   int i;
};


class runner1 : public runnable {
public:
   runner1(queue<item, synchronized> *q) {
      this->q = q;
   }
   
   virtual void run() {
      int n = 0;
      for (n = 0; n < 10; n++) {
         this->q->enqueue(*(new item(n)));
         /* std::cout << n << "++\n"; */
         thread::sleep_usec(10000);
      }
   }
protected:
   queue<item, synchronized> *q;
};
      

class runner2 : public runnable {
public:
   runner2(queue<item, synchronized> *q) {
      this->q = q;
   }
   
   virtual void run() {
      int n = 0;
      for (n = 10; n > 0; n--) {
         this->q->enqueue(*(new item(n)));
         /* std::cout << n << "--\n"; */
         thread::sleep_nsec(10);
      }
   }
protected:
   queue<item, synchronized> *q;
};

int main() {
   queue<item, synchronized> *q = new queue<item, synchronized>(compare_functors<item>::f_compare());
   runnable *r1 = new runner1(q);
   runnable *r2 = new runner2(q);
   thread *t1 = new thread((*r1));
   thread *t2 = new thread((*r2));
   t1->start();
   t2->start();
   t1->join();
   t2->join();
   std::cout << "dequeueing...";
   while (!(q->is_empty())) {
      q->dequeue().print();
   }
   std::cout << "\n";
   delete t1;
   delete t2;
   delete r1;
   delete r2;
   delete q;
   return 0;
}
