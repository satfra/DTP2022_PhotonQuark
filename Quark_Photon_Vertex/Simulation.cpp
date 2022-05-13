#include <iostream>
#include <LegendrePolynomials.hh>
#include <QuadratureIntegral.hh>
#include <Utils.hh>
#include <LinearInterpolate.hh>
#include <iteration.hh>
#include <numeric>

int main(int argc, char *argv[]) 
{
  constexpr unsigned order_k = 256;
  constexpr unsigned order_q = 1;
  constexpr unsigned order_z = 2;
  static_assert(order_z % 2 == 0);

  std::vector<double> k_grid(order_k);
  std::vector<double> q_grid(order_q);

  std::iota(k_grid.begin(), k_grid.end(), 0);
  std::iota(q_grid.begin(), q_grid.end(), 0);

  constexpr double L = 5e+1;
  std::cout << "L=" << L << " order_k=" << order_k <<"\n";
  constexpr double k_sq_min = 1e-12;
  constexpr double k_sq_max = L*L;
  k_grid = linearMapTo(k_grid, 0., double(k_grid.size()-1), k_sq_min, k_sq_max);
  q_grid = linearMapTo(q_grid, 0., double(q_grid.size()-1), 1e-3, 3.);

  q_grid[0] = 1.;

  LegendrePolynomial<order_z> lp;
  const std::vector<double> z_grid = lp.zeroes();
  const std::vector<double> y_grid = z_grid;

  iterate_a_and_b(q_grid, z_grid, k_grid, y_grid);

  return 0;
}
