/*
 * Test throwing an exception.
 */
#include "lang/exceptions/exception.hh"
#include <iostream>

using lang::exceptions::exception;

int main() {
   try {
      throw (exception("an exception"));
   } catch (exception& e) {
      std::cout << e.what() << "\n";
   }
   return 0;
}
