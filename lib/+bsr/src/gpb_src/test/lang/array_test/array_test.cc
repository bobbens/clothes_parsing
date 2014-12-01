/*
 * Test array operations.
 */
#include <iostream>
#include "lang/array.hh"

using lang::array;

array<double> get_array() {
   array<double> a(20);
   a(0) = 0;
   a(1) = 0.166035;
   a(2) = 0.0381288;
   a(3) = 0.00558394;
   a(4) = 0.191116;
   a(5) = 0.0099273;
   a(6) = 0.198789;
   a(7) = 0.23788;
   a(8) = 0.1603;
   a(9) = 0.264449;
   a(10) = 0.297406;
   a(11) = 0.0491625;
   a(12) = 0.155613;
   a(13) = 0.346112;
   a(14) = 0.119747;
   a(15) = 0.340048;
   a(16) = 0.13701;
   a(17) = 0.314217;
   a(18) = 0;
   a(19) = 0;
   return a;
}

void test_sort() {
   array<double> a = get_array();
   std::cout << "unsorted = ";
   for (unsigned long n = 0; n < a.size(); n++)
      std::cout << a(n) << " ";
   std::cout << "\n";
   a.sort();
   std::cout << "sorted   = ";
   for (unsigned long n = 0; n < a.size(); n++)
      std::cout << a(n) << " ";
   std::cout << "\n";
}

void test_sort_idx() {
   array<double> a = get_array();
   std::cout << "unsorted = ";
   for (unsigned long n = 0; n < a.size(); n++)
      std::cout << a(n) << " ";
   std::cout << "\n";
   array<unsigned long> idx = a.sort_idx();
   std::cout << "sorted   = ";
   for (unsigned long n = 0; n < a.size(); n++)
      std::cout << a(n) << " ";
   std::cout << "\n";
   std::cout << "idx = ";
   for (unsigned long n = 0; n < idx.size(); n++)
      std::cout << idx(n) << " ";
   std::cout << "\n";
}

void test_unique() {
   array<double> a = get_array();
   std::cout << "array    = ";
   for (unsigned long n = 0; n < a.size(); n++)
      std::cout << a(n) << " ";
   std::cout << "\n";
   a.unique();
   std::cout << "unique   = ";
   for (unsigned long n = 0; n < a.size(); n++)
      std::cout << a(n) << " ";
   std::cout << "\n";
}

void test_unique_idx() {
   array<double> a = get_array();
   std::cout << "array    = ";
   for (unsigned long n = 0; n < a.size(); n++)
      std::cout << a(n) << " ";
   std::cout << "\n";
   array<unsigned long> idx = a.unique_idx();
   std::cout << "unique   = ";
   for (unsigned long n = 0; n < a.size(); n++)
      std::cout << a(n) << " ";
   std::cout << "\n";
   std::cout << "idx = ";
   for (unsigned long n = 0; n < idx.size(); n++)
      std::cout << idx(n) << " ";
   std::cout << "\n";
}

int main() {
   test_sort();
   test_sort_idx();
   test_unique();
   test_unique_idx();
   return 0;
}
