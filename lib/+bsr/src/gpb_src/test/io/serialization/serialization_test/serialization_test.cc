/*
 * Serialization test.
 */
#include "collections/array_list.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "io/serialization/serial_file_istream.hh"
#include "io/serialization/serial_file_ostream.hh"
#include "io/streams/cout.hh"
#include "lang/exceptions/exception.hh"
#include "lang/pointers/auto_ptr.hh"
#include "lang/string.hh"
#include "math/complex.hh"
#include "math/matrices/cmatrix.hh"

using collections::array_list;
using collections::list;
using collections::pointers::auto_collection;
using io::serialization::serial_file_istream;
using io::serialization::serial_file_ostream;
using io::streams::cout;
using lang::exceptions::exception;
using lang::pointers::auto_ptr;
using lang::string;
using math::complex;
using math::matrices::cmatrix;

int main() {
   try {
      {
         cmatrix<> m(5,5);
         for (unsigned long n = 0; n < m.size(); n++)
            m[n] = complex<>(n, n);
         complex<> z(256,78);
         serial_file_ostream os("results/serialization_test/serial_file");
         z.serialize(os);
         string<>("hi").serialize(os);
         m.serialize(os);
         list< cmatrix<> > lst;
         cmatrix<> mt = transpose(m);
         lst.add(m);
         lst.add(mt);
         lst.serialize(os);
         array_list< cmatrix<> > a(lst);
         a.serialize(os);
      }
      {
         serial_file_istream is("results/serialization_test/serial_file");
         auto_ptr< complex<> > z = complex<>::deserialize(is);
         cout << *z << "\n";
         auto_ptr< string<> > s = string<>::deserialize(is);
         cout << *s << "\n";
         auto_ptr< cmatrix<> > m = cmatrix<>::deserialize(is);
         cout << *m << "\n";
         auto_collection< cmatrix<>, list< cmatrix<> > > lst =
            list< cmatrix<> >::deserialize(is);
         cout << lst->head() << "\n";
         cout << lst->tail() << "\n";
         auto_collection< cmatrix<>, array_list< cmatrix<> > > a =
            array_list< cmatrix<> >::deserialize(is);
         cout << a->head() << "\n";
         cout << a->tail() << "\n";
      }
   } catch (exception& e) {
      cout << e << "\n";
   }
   return 0;
}
