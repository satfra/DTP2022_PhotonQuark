#include <iostream>
#include <LegendrePolynomials.hh>
#include <QuadratureIntegral.hh>
#include <Utils.hh>
#include <LinearInterpolate.hh>
#include <iteration.hh>
#include <numeric>
#include <parameters.hh>

int main(int argc, char *argv[]) 
{
  static_assert(parameters::numerical::z_steps % 2 == 0);

  std::vector<double> k_grid(parameters::numerical::k_steps);
  std::vector<double> q_grid(parameters::numerical::q_steps);

  std::iota(k_grid.begin(), k_grid.end(), 0);
  std::iota(q_grid.begin(), q_grid.end(), 0);

  std::cout << "L=" << parameters::physical::lambda_UV << " order_k=" << parameters::numerical::k_steps <<"\n";
  k_grid = linearMapTo(k_grid, 0., double(k_grid.size()-1), std::log(parameters::physical::lambda_IR),
                       std::log(parameters::physical::lambda_UV));
  q_grid = linearMapTo(q_grid, 0., double(q_grid.size()-1), 1e-3, 3.);

  LegendrePolynomial<parameters::numerical::z_steps> lp_z;
  LegendrePolynomial<parameters::numerical::y_steps> lp_y;
  const std::vector<double> z_grid = lp_z.zeroes();
  const std::vector<double> y_grid = lp_y.zeroes();

  iterate_a_and_b(q_grid, z_grid, k_grid, y_grid);

  return 0;
}
