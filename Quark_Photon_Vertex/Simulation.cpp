#include <iostream>
#include <LegendrePolynomials.hh>
#include <QuadratureIntegral.hh>
#include <Utils.hh>
#include <LinearInterpolate.hh>
#include <iteration.hh>
#include <numeric>

int main(int argc, char *argv[]) {
    // constexpr unsigned order_Q = 10;
    // constexpr unsigned order_k = 10;
/*
  qIntegral2d<LegendrePolynomial<order_k>, LegendrePolynomial<order_Q>> qint;
  auto f = [](const double& k, const double& Q){ return 2.; };
  std::cout << "We get " << qint(f, 0., 1., 0., 2.) << "\n";*/

    constexpr unsigned int order = 4;
    static_assert(order % 2 == 0);
    std::vector<double> k_grid(order);
    std::vector<double> q_grid(order);

    std::iota(k_grid.begin(), k_grid.end(), 0);
    std::iota(q_grid.begin(), q_grid.end(), 0);

    const double k_min = 1e-5;
    const double k_max = 1e6;
    k_grid = linearMapTo(k_grid, 0., double(k_grid.size()-1), k_min, k_max);
    q_grid = linearMapTo(q_grid, 0., double(q_grid.size()-1), 1e-3, 3.);

    for (unsigned i = 0; i < q_grid.size(); ++i)
      std::cout << q_grid[i] << " q ";
    std::cout << "\n";

    LegendrePolynomial<order> lp;
    std::vector<double> z_grid = lp.zeroes();
    for (unsigned i = 0; i < z_grid.size(); ++i)
      std::cout << z_grid[i] << " z ";
    std::cout << "\n";
    const std::vector<double>& y_grid = z_grid;
    iterate_a_and_b(q_grid, z_grid, k_grid, y_grid);

    return 0;
}
