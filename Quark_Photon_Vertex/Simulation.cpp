#include <iostream>
#include <numeric>

#include "Utils.hh"
#include "LegendrePolynomials.hh"
#include "quark_model_functions.hh"
#include "iteration.hh"
#include "parameters.hh"

int main(int argc, char *argv[]) 
{
  // get flags from shell
  std::string flags = argc > 1 ? argv[1] : "";

  const bool debug = flags.find('v') < flags.length() ? true : false;
  if(debug) std::cout << "Showing debug output.\n";

  const bool use_quark_DSE = flags.find('d') < flags.length() ? true : false;
  if(use_quark_DSE) std::cout << "Using the quark DSE.\n";

  const bool use_PauliVillars = flags.find('p') < flags.length() ? true : false;
  if(use_PauliVillars) std::cout << "Using Pauli-Villars regularisation.\n";

  // avoid z == 0 in a grid, which would lead to division by zero.
  static_assert(parameters::numerical::z_steps % 2 == 0);

  // create the k_grid, fill it and transform it to the correct range
  std::vector<double> k_grid(parameters::numerical::k_steps);
  std::iota(k_grid.begin(), k_grid.end(), 0);
  k_grid = linearMapTo(k_grid, 0., double(k_grid.size()-1),
                       std::log(parameters::physical::lambda_IR),
                       std::log(parameters::physical::lambda_UV));

  // create the q_grid, fill it and transform it to the correct range
  std::vector<double> q_grid(parameters::numerical::q_steps);
  std::vector<double> q_gridtemp(parameters::numerical::q_steps);
  std::iota(q_grid.begin(), q_grid.end(), 0);

  q_gridtemp = linearMapTo(q_gridtemp, 0., double(q_grid.size()-1), std::log(parameters::numerical::min_q_sq),
                                std::log(parameters::numerical::max_q_sq));
  for (unsigned int i = 0; i < q_grid.size(); ++i) {
    q_grid[i] = std::exp(q_gridtemp[i]);
  }

  // We use the zeroes of LegendrePolynomials for the z and y grids
  LegendrePolynomial<parameters::numerical::z_steps> lp_z;
  LegendrePolynomial<parameters::numerical::y_steps> lp_y;
  const std::vector<double> z_grid = lp_z.zeroes();
  const std::vector<double> y_grid = lp_y.zeroes();

  // Start the program
  if(use_quark_DSE)
    iterate_a_and_b<quark_DSE>(q_grid, z_grid, k_grid, y_grid, use_PauliVillars, debug);
  else
    iterate_a_and_b<quark_model>(q_grid, z_grid, k_grid, y_grid, use_PauliVillars, debug);

  return 0;
}
