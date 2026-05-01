#pragma once

#include <cmath>

#include "Utils.hh"

// Reduced grid sizes for performance benchmarking. Selected via the CMake
// option -DBENCHMARK_PROFILE=ON. Production parameters live in parameters.hh.
//
// Sizes are chosen so each flag combo runs in ~2-3 min (fast enough to iterate
// per-tier) while staying large enough that timing variance is well below the
// 5-15% tier-level gains we expect to measure. Cost scaling vs production:
//   K-kernel ∝ n_structs² · k² · z² · y → factor (128/64)² · (32/24)² · (32/24)
//                                       ≈ 9.5× faster per q-iteration
//   q-loop  → 32/8 = 4× fewer q-points
//   total   ≈ 38× faster wall-clock.

namespace parameters {
namespace numerical {
constexpr unsigned k_steps = 64;
constexpr unsigned z_steps = 24;
constexpr unsigned y_steps = 24;
constexpr unsigned q_steps = 8;

constexpr double min_q_sq = 1e-5;
constexpr double max_q_sq = 1;

constexpr double target_acc = 1e-5;
constexpr unsigned max_steps = 40;

constexpr unsigned quark_dse_steps_q = 400;
constexpr unsigned quark_dse_steps_z = 128;
constexpr unsigned quark_dse_max_steps = 80;
constexpr double quark_dse_acc = 1e-8;

constexpr unsigned int n_structs = 12;
constexpr double int_factors = 0.5 / powr<4>(2. * M_PI);
} // namespace numerical

namespace physical {
constexpr double lambda_UV = 1e6;
constexpr double lambda_IR = 1e-6;

constexpr double eta_mt = 1.8;
constexpr double lambda_mt = 0.72;

constexpr double lambda_qcd = 0.234;
constexpr double lambda_0 = 1.0;
constexpr double gamma_m = 0.48;

constexpr double lambda_pv = 200.0;

constexpr double m_c = 0.0037;
constexpr double mu = 19.0;
constexpr double quark_a0 = 1.0;
} // namespace physical
} // namespace parameters
