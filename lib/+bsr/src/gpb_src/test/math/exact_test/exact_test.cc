/*
 * Test exact arithmetic.
 */
#include "io/streams/cout.hh"
#include "math/exact.hh"
#include "math/math.hh"

using io::streams::cout;
using math::exact;

template <typename T>
void test_sum(const T& small, const T& large) {
   T sum = small + large;
   cout << "sum = " << sum << "\n";
   cout << "sum - small = " << (sum - small) << "\n";
   cout << "sum - large = " << (sum - large) << "\n";
}

template <typename T>
void test_diff(const T& small, const T& large) {
   T diff = large - small;
   cout << "diff = " << diff << "\n";
   cout << "diff + small = " << (diff + small) << "\n";
   cout << "large - diff = " << (large - diff) << "\n";
}

template <typename T>
void test_prod(const T& t0, const T& t1) {
   T prod = t0 * t1;
   cout << "prod = " << prod << "\n";
}

void test_exact() {
   exact<> a(0.25);
   exact<> b(1244);
   exact<> pi(M_PIl);
   exact<> pi_2 = 0.5 * pi;
   exact<> q(2.234234e+156);
   exact<> r(0.234234e-23);
   exact<> x = pi_2*a + pi*b + q + 5.34 + a*b + pi*pi_2 + r;
   cout << "x = " << x << " [" << x.size() << " term(s)]" << "\n";
   exact<> y = x;
   cout << "compressing x\n";
   x.compress();
   cout << "x = " << x << " [" << x.size() << " term(s)]" << "\n";
   cout << "compresson " << ((x == y) ? "passed" : "failed") << "\n";
   exact<> ab = x - pi*pi_2 - pi_2*a - 5.34 - pi*b - q - r;
   cout << "a*b = " << ab << " [" << ab.size() << " term(s)]" << "\n";
   cout << "compressing ab\n";
   ab.compress();
   cout << "a*b = " << ab << " [" << ab.size() << " term(s)]" << "\n";
   cout << "a*b == 311 " << ((ab == 311.0) ? "(true)" : "(false)") << "\n";
}

int main() {
   cout << "epsilon (relative error bound) = " << exact<>::epsilon() << "\n\n";
   cout << "--- double precision arithmetic ---\n";
   double small = 1.203984e-204;
   double large = 1.203984e+204;
   cout << "small = " << small << "\n";
   cout << "large = " << large << "\n";
   test_sum(small, large);
   test_diff(small, large);
   test_prod(small, large);
   cout << "\n--- exact arithmetic ---\n";
   exact<> e_small(small);
   exact<> e_large(large);
   cout << "small = " << e_small << "\n";
   cout << "large = " << e_large << "\n";
   test_sum(e_small, e_large);
   test_diff(e_small, e_large);
   test_prod(e_small, e_large);
   cout << "\n--- more exact arithmetic ---\n";
   test_exact();
   return 0;
}
