#include <iostream>
#include <numeric>

#include "Utils.hh"
#include "LegendrePolynomials.hh"
#include "ChebyshevPolynomial2.hh"
#include "quark_model_functions.hh"
#include "iteration.hh"
#include "hvp.hh"
#include "flavor.hh"
#include "parameters.hh"

int main(int argc, char *argv[])
{
  // get flags from shell
  std::string flags = argc > 1 ? argv[1] : "";

  const bool debug = flags.find('v') != std::string::npos;
  if(debug) std::cout << "Showing debug output.\n";

  const bool use_quark_DSE = flags.find('d') != std::string::npos;
  if(use_quark_DSE) std::cout << "Using the quark DSE.\n";

  const bool use_PauliVillars = flags.find('p') != std::string::npos;
  if(use_PauliVillars) std::cout << "Using Pauli-Villars regularisation.\n";

  const bool calculate_hvp = flags.find('h') != std::string::npos;
  if(calculate_hvp) std::cout << "Calculating the hadronic vacuum polarisation after the QPV.\n";

  // Quark flavour: first of l/s/c/b found in the flag string (default light).
  char flavor_char = 'l';
  for (char c : {'l', 's', 'c', 'b'})
    if (flags.find(c) != std::string::npos) { flavor_char = c; break; }
  const FlavorParams fp = flavor_from_flag(flavor_char);
  std::cout << "Quark flavour: " << fp.name
            << " (m_c = " << fp.m_c << " GeV, eta = " << fp.eta_mt
            << ", Lambda = " << fp.lambda_mt << " GeV)\n";

  // The analytic quark_model is a hardcoded light-quark fit and ignores m_c,
  // so any non-light flavour is meaningless without the quark DSE (-d).
  if (fp.name != "light" && !use_quark_DSE) {
    std::cerr << "\nError: flavour '" << fp.name << "' requires the quark DSE. "
              << "Add the 'd' flag (the analytic quark model ignores the quark "
              << "mass and is only valid for the light quark).\n";
    return 1;
  }

  // One HDF5 file per run, tagged by flavour + mode.
  const std::string run_tag = fp.name
      + (use_quark_DSE ? "_dse" : "")
      + (use_PauliVillars ? "_pv" : "");
  const std::string h5_path = "output_" + run_tag + ".h5";

  // avoid z == 0 in a grid, which would lead to division by zero.
  static_assert(parameters::numerical::z_steps % 2 == 0);

  // create the k_grid in log(k²), optionally biased toward the IR: the
  // uniform parameter u = i/(N−1) ∈ [0,1] is remapped through the cubic
  // f(u) = α·u + (1−α)·u³ with α = k_grid_ir_bias (see parameters.hh).
  std::vector<double> k_grid(parameters::numerical::k_steps);
  {
    using parameters::numerical::k_grid_ir_bias;
    const double log_IR = std::log(parameters::physical::lambda_IR);
    const double log_UV = std::log(parameters::physical::lambda_UV);
    const double L = log_UV - log_IR;
    const double n_minus_1 = double(k_grid.size() - 1);
    for (std::size_t i = 0; i < k_grid.size(); ++i) {
      const double u = double(i) / n_minus_1;
      const double f = k_grid_ir_bias * u + (1.0 - k_grid_ir_bias) * u * u * u;
      k_grid[i] = log_IR + f * L;
    }
  }

  // create the q_grid, fill it and transform it to the correct range
  std::vector<double> q_grid(parameters::numerical::q_steps);
  std::iota(q_grid.begin(), q_grid.end(), 0);
  q_grid = linearMapTo(q_grid, 0., double(q_grid.size()-1),
                       std::log(parameters::numerical::min_q_sq),
                       std::log(parameters::numerical::max_q_sq));
  std::transform(q_grid.begin(), q_grid.end(), q_grid.begin(),
                 [](double x) { return std::exp(x); });

  // y is a uniform-weight angular integration → Legendre.
  // z carries a √(1−z²) Jacobian from the 4D measure → Chebyshev type 2,
  // which absorbs that weight by construction.
  ChebyshevPolynomial2<parameters::numerical::z_steps> cp_z;
  LegendrePolynomial<parameters::numerical::y_steps> lp_y;
  const std::vector<double> z_grid = cp_z.zeroes();
  const std::vector<double> y_grid = lp_y.zeroes();

  // Start the program
  HvpBInput hvp_in = use_quark_DSE
    ? iterate_a_and_b<quark_DSE>(q_grid, z_grid, k_grid, y_grid, use_PauliVillars, debug, fp, use_quark_DSE, h5_path)
    : iterate_a_and_b<quark_model>(q_grid, z_grid, k_grid, y_grid, use_PauliVillars, debug, fp, use_quark_DSE, h5_path);

  if(calculate_hvp)
    hvp::hvp_driver(hvp_in, q_grid, k_grid, z_grid, h5_path);

  return 0;
}
