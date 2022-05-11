#include <iostream>
#include <LegendrePolynomials.hh>
#include <QuadratureIntegral.hh>
#include <Utils.hh>
#include <LinearInterpolate.hh>
#include <iteration.hh>

int main(int argc, char *argv[]) {
    // constexpr unsigned order_Q = 10;
    // constexpr unsigned order_k = 10;
/*
  qIntegral2d<LegendrePolynomial<order_k>, LegendrePolynomial<order_Q>> qint;
  auto f = [](const double& k, const double& Q){ return 2.; };
  std::cout << "We get " << qint(f, 0., 1., 0., 2.) << "\n";*/

    constexpr unsigned int order = 2;
    LegendrePolynomial<order> lp;
    const auto grid1 = lp.zeroes();
    iterate_a_and_b(grid1, grid1, grid1, grid1);

    return 0;
}
