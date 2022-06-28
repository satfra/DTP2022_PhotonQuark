#pragma once

#include <cmath>

#include "Utils.hh"

namespace parameters
{
  namespace numerical
  {
    // The number of steps in the k/k' grid
    constexpr unsigned k_steps = 64;
    // The number of steps in the z/z' grid
    constexpr unsigned z_steps = 16;
    // The number of steps in the y grid
    constexpr unsigned y_steps = 16;
    // The number of steps in the Q grid
    constexpr unsigned q_steps = 2;

    constexpr double min_q_sq = 1e-4;
    constexpr double max_q_sq = 1;

    // Target accuracy for the iteration
    constexpr double target_acc = 1e-5;
    // Maximum number of iteration steps
    constexpr unsigned max_steps = 100;

    // Grid for the quark propagator dse
    constexpr unsigned quark_dse_steps_q = 512;
    constexpr unsigned quark_dse_steps_z = 64;
    constexpr unsigned quark_dse_max_steps = 200;
    constexpr double quark_dse_acc = 1e-8;

    // number of tensor structures, is ALWAYS fixed to 12, don't change
    constexpr unsigned int n_structs = 12;
    constexpr unsigned int n_structs_L = 4;
    constexpr unsigned int n_structs_T = 8;
    // integration factor, don't change
    constexpr double int_factors = 0.5 / powr<4>(2. * M_PI);
  }

  namespace physical
  {
     // The IR and UV cutoff for k^2 in the BSE
    constexpr double lambda_UV = 1e6;
    constexpr double lambda_IR = 1e-4;

    // The IR and UV cutoff for k^2 in the DSE
    constexpr double lambda_UV_DSE = 1e6;
    constexpr double lambda_IR_DSE = 1e-6;
    
    // The IR parameters for Maris-Tandy
    constexpr double eta_mt = 1.8;
    constexpr double lambda_mt = 0.72;

    // The UV parameters for Maris-Tandy
    constexpr double lambda_qcd = 0.234;
    constexpr double lambda_t = 1.0;
    constexpr double gamma_m = 0.48;

    // Scale for Pauli-Villars
    constexpr double lambda_pv = 200.0;

    // Parameters for the quark dse
    constexpr double m_c = 0.0038;
    constexpr double mu = 19.0;
    constexpr double quark_a0 = 1.0;
  }
}
