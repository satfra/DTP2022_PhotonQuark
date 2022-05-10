#include <iostream>
#include <LegendrePolynomials.hh>
#include <QuadratureIntegral.hh>
#include <Utils.hh>
#include <LinearInterpolate.hh>

int main(int argc, char* argv[])
{
  constexpr unsigned order_Q = 10;
  constexpr unsigned order_k = 10;
/*
  qIntegral2d<LegendrePolynomial<order_k>, LegendrePolynomial<order_Q>> qint;
  auto f = [](const double& k, const double& Q){ return 2.; };
  std::cout << "We get " << qint(f, 0., 1., 0., 2.) << "\n";*/

  std::vector<double> x{1,2,3,4,5};
  std::vector<double> f{1,9,3,4,5};
  lInterpolator lpol(x, f);
  std::cout << "Interp is " <<lpol(1.5) << "\n";

  return 0;
}
