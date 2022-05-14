#pragma once

namespace parameters {
    namespace numerical {
        // The number of steps in the k/k' grid
        constexpr unsigned int k_steps = 64;
        // The number of steps in the z/z' grid
        constexpr unsigned int z_steps = 16;
        // The number of steps in the y grid
        constexpr unsigned int y_steps = 12;
        // The number of steps in the Q grid
        constexpr unsigned int q_steps = 2;

        // Target accuracy for the iteration
        constexpr double target_acc = 1e-3;
        // Maximum number of iteration steps
        constexpr unsigned max_steps = 10;
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

        // The UV parameters for Maris-Tandy
        constexpr double lambda_qcd = 0.234;
        constexpr double lambda_0 = 1.0;
        constexpr double gamma_m = 0.48;
    }
}
