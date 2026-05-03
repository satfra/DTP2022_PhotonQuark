#pragma once

#ifdef BENCHMARK_PROFILE
#include "parameters_bench.hh"
#else

#include <cmath>

#include "Utils.hh"

namespace parameters {
namespace numerical {
// The number of steps in the k/k' grid
constexpr unsigned k_steps = 128;
// The number of steps in the z/z' grid
constexpr unsigned z_steps = 32;
// The number of steps in the y grid
constexpr unsigned y_steps = 32;
// The number of steps in the Q grid
constexpr unsigned q_steps = 32;

constexpr double min_q_sq = 1e-5;
constexpr double max_q_sq = 1;

// Target accuracy for the iteration
constexpr double target_acc = 1e-5;
// Maximum number of iteration steps
constexpr unsigned max_steps = 100;

// IR bias for the log(k²) grid built in Simulation.cpp. Maps u ∈ [0,1]
// through f(u) = α·u + (1−α)·u³ with α = k_grid_ir_bias. α = 1 reproduces
// the uniform-log grid; α < 1 packs more points near lambda_IR. α = 0.3
// gives ≈3.3× IR density and ≈2.06× as many points in the first three
// decades above lambda_IR (where the WTI2 residual lives), at the cost of
// a 1.32× wider UV step than the uniform grid.
constexpr double k_grid_ir_bias = 0.3;

// Grid for the quark propagator dse
constexpr unsigned quark_dse_steps_q = 1000;
constexpr unsigned quark_dse_steps_z = 256;
constexpr unsigned quark_dse_max_steps = 200;
constexpr double quark_dse_acc = 1e-8;

// number of tensor structures, is ALWAYS fixed to 12, don't change
constexpr unsigned int n_structs = 12;
// integration factor, don't change
constexpr double int_factors = 0.5 / powr<4>(2. * M_PI);
} // namespace numerical

namespace physical {
// The UV cutoff for k^2
constexpr double lambda_UV = 1e6;
// The IR cutoff for k^2
constexpr double lambda_IR = 1e-6;

// The IR parameters for Maris-Tandy
constexpr double eta_mt = 1.8;
constexpr double lambda_mt = 0.72;

// The UV parameters for Maris-Tandy
constexpr double lambda_qcd = 0.234;
constexpr double lambda_0 = 1.0;
constexpr double gamma_m = 0.48;

// Scale for Pauli-Villars
constexpr double lambda_pv = 200.0;

// Parameters for the quark dse
constexpr double m_c = 0.0037;
constexpr double mu = 19.0;
constexpr double quark_a0 = 1.0;
} // namespace physical
} // namespace parameters

#endif // BENCHMARK_PROFILE
