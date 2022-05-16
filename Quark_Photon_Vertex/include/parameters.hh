#pragma once

namespace parameters {
    namespace numerical {
        // number of tensor structures, is ALWAYS fixed to 12
        constexpr unsigned int n_structs = 12;
        // The number of steps in the k/k' grid
        constexpr unsigned int k_steps = 32;
        // The number of steps in the z/z' grid
        constexpr unsigned int z_steps = 16;
        // The number of steps in the y grid
        constexpr unsigned int y_steps = 12;
        // The number of steps in the Q grid
        constexpr unsigned int q_steps = 2;

        // Target accuracy for the iteration
        constexpr double target_acc = 1e-3;
        // Maximum number of iteration steps
        constexpr unsigned max_steps = 30;
        
        // integration factor, should be fixed
        constexpr double int_factors = 0.5 / powr<4>(2. * M_PI);

        // Grid for the quark propagator dse
        constexpr unsigned int quark_dse_steps_q = 600;
        constexpr unsigned int quark_dse_steps_z = 128;
        constexpr unsigned int quark_dse_max_steps = 100;
        constexpr double quark_dse_acc = 1e-7;
    }

    namespace physical {
        // The UV cutoff for k^2
        constexpr double lambda_UV = 1e6;
        // The IR cutoff for k^2
        constexpr double lambda_IR = 1e-6;

        // The value for z_2
        constexpr double z_2 = 0.97;

        // The parameters for Maris-Tandy
        constexpr double eta_mt = 1.8;
        constexpr double lambda_mt = 0.72;
        constexpr bool pauliVillars = true;

        // The UV parameters for Maris-Tandy
        constexpr double lambda_qcd = 0.234;
        constexpr double lambda_0 = 1.0;
        constexpr double gamma_m = 0.48;
    }
}
