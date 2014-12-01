/*
 * Test fractions.
 */
#include "io/streams/cout.hh"
#include "math/fraction.hh"

using io::streams::cout;
using math::fraction;

int main() {
   fraction<int> z(0,10);
   fraction<int> a(1,12);
   fraction<int> b(3,12);
   fraction<int> c(7,23);
   fraction<double> f(4,12);
   cout << z << "\n";
   cout << z.reduce() << "\n";
   cout << b.reduce() << "\n";
   cout << f << "\n";
   cout << f.reduce() << "\n";
   cout << "a = " << a << "\n";
   cout << "b = " << b << "\n";
   cout << "c = " << c << "\n";
   cout << "a + b = " << a + b << "\n";
   cout << "a + c = " << a + c << "\n";
   cout << "a - c = " << a - c << "\n";
   cout << "a * c = " << a * c << "\n";
   cout << "a / c = " << a / c << "\n";
   return 0;
}
