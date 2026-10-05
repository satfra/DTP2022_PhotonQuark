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
// The number of steps in the y grid (gluon angle). WTIs are converged from
// y ≈ 12 on; 32 leaves margin. k_steps = 128 is what the bottom quark needs.
constexpr unsigned y_steps = 32;
// The number of steps in the Q grid
constexpr unsigned q_steps = 48;

// Π(Q²) − Π(0) is linear in Q² below ~1e-2 GeV²; points far below that only
// anchor the Π(0) fit in hvp.hh, and below ~1e-4 they are noise-dominated.
constexpr double min_q_sq = 1e-4;
constexpr double max_q_sq = 4;
// Upper end of the low-p² fit that extracts Π(0) for the HVP (hvp.hh);
// needs ≥ 6 q-grid points below it.
constexpr double hvp_fit_max_p_sq = 1e-2;

// Target accuracy for the iteration
constexpr double target_acc = 1e-7;
// Maximum number of iteration steps
constexpr unsigned max_steps = 100;

// IR bias for the log(k²) grid built in Simulation.cpp. Maps u ∈ [0,1]
// through f(u) = α·u + (1−α)·u³ with α = k_grid_ir_bias. α = 1 reproduces
// the uniform-log grid; α < 1 packs more points near lambda_IR. Scans showed
// no WTI gain from α < 1, so the uniform grid is the default.
constexpr double k_grid_ir_bias = 1.0;

// Grid for the quark propagator dse
constexpr unsigned quark_dse_steps_q = 1024;
constexpr unsigned quark_dse_steps_z = 256;
constexpr unsigned quark_dse_max_steps = 500;
constexpr double quark_dse_acc = 1e-9;

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
