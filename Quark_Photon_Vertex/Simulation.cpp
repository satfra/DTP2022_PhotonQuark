#include <iostream>
#include <LegendrePolynomials.hh>
#include <QuadratureIntegral.hh>
#include <Utils.hh>
#include <LinearInterpolate.hh>

  template<unsigned order>
void checkForOrder()
{
  auto f1 = [](const double& x){ return std::sqrt(1. - powr<2>(x)); };
  constexpr double f1_a = -1.;
  constexpr double f1_b = 1.;
  constexpr double f1_should = M_PI / 2.;

  auto f2 = [](const double& x){ return powr<2>(std::log(1. + x)); };
  constexpr double f2_a = 0.;
  constexpr double f2_b = 1.;
  constexpr double f2_should = 2. * powr<2>(std::log(2.)) - 4. * std::log(2.) + 2.;

  auto f3 = [](const double& x){ return std::log(1. - x) / x; };
  constexpr double eps = 1e-10;
  constexpr double f3_a = 0. + eps;
  constexpr double f3_b = 1. - eps;
  constexpr double f3_should = - powr<2>(M_PI) / 6.;

  auto f4 = [](const double& x){ return powr<2>(std::log(1. - x)); };
  constexpr double f4_a = 0.;
  constexpr double f4_b = 1.;
  constexpr double f4_should = 2.;

  qIntegral<LegendrePolynomial<order>> qint;

  auto compareQIntegal = [&](auto fun, auto a, auto b, const auto should) {
    auto res = qint(fun, a, b);
    std::cout << "should: " << should << ", is: " << res << ", diff: " << std::abs(should-res)/should << "\n";
  };
  compareQIntegal(f1, f1_a, f1_b, f1_should);
  compareQIntegal(f2, f2_a, f2_b, f2_should);
  compareQIntegal(f3, f3_a, f3_b, f3_should);
  compareQIntegal(f4, f4_a, f4_b, f4_should);
}

template<unsigned order = 0>
void cfo(unsigned o)
{
  constexpr unsigned MAXORDER = 150;
  if (order == o)
    checkForOrder<order>();
  else if constexpr (order < MAXORDER)
    cfo<order+1>(o);
  else
    throw std::runtime_error("nope");
}

int main(int argc, char* argv[])
{
  cfo(100);

  auto fun = [](const auto& x){ return 0.2/(x + 0.5) + 0.002;};
  constexpr unsigned N = 10;

  std::vector<double> x(N);
  std::vector<double> z(N);
  std::vector<double> f(N);

  for(unsigned i = 0; i < N; ++i)
  {
    z[i] = -1. + 2.*i / double(N);
    f[i] = fun(z[i]);
    std::cout << z[i] << " - " << f[i] << "\n";
  }

  auto z1(x); 
  std::for_each(z1.begin(), z1.end(), [](auto& v){ v = powr<1>((v-1)/(v+1)); });

  return 0;
}
