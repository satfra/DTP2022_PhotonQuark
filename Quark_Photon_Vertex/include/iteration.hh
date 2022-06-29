#pragma once

#include <complex>

// headers from /include
#include "types.hh"
#include "parameters.hh"
#include "Kernels_K.hh"
#include "Kernels_G.hh"
#include "QuadratureIntegral.hh"
#include "LegendrePolynomials.hh"
#include "Utils.hh"
#include "fileIO.hh"
#include "maris_tandy.hh"
#include "basistransform.hh"
#include "WTI.hh"
#include "kernel_calculation.hh"
#include <chrono>
// #include "quark_dressings.hh"

using namespace std;
using namespace parameters::numerical;

using Integrator1d = qIntegral<LegendrePolynomial<parameters::numerical::y_steps>>;
using DiscrIntegrator1d = yIntegral<LegendrePolynomial<parameters::numerical::y_steps>>;
using Integrator2d = qIntegral2d<LegendrePolynomial<parameters::numerical::k_steps>, LegendrePolynomial<parameters::numerical::z_steps>>;
using DiscrIntegrator2d = gaulegIntegral2d<LegendrePolynomial<parameters::numerical::k_steps>, LegendrePolynomial<parameters::numerical::z_steps>>;


double a0 (const unsigned& i)
{
  if(i == 0)
    return sqrt(2.);
  else if(i == 6 || i == 9)
    return 1.;
  return 0.;
}

double jakobian_log_k_prime (const double &z_prime, const double &log_k_prime_sq)
{
    return parameters::numerical::int_factors * 2. * M_PI * sqrt(1. - powr<2>(z_prime)) * std::exp(2 * log_k_prime_sq);
}

void a_initialize(
    dressing* a, // BSE dressings a
    const double z2 // quark renormalization constant
    )
{
  using namespace parameters::numerical;
  for (unsigned i = 0; i < n_structs; ++i)
  {
    #pragma omp parallel for collapse(2)
    for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
      for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx)
        (*a)[i][k_idx][z_idx] = z2 * a0(i);
  }
}


mat_cmplx average_array_z0(const dressing* a, const unsigned z_0)
{
  using namespace parameters::numerical;

  mat_cmplx a_z0((*a).size(), vec_cmplx((*a)[0].size(), 0.0));

  for (unsigned i = 0; i < (*a).size(); ++i)
    for (unsigned k_idx = 0; k_idx < (*a)[i].size(); ++k_idx)
      a_z0[i][k_idx] = 0.5 * ((*a)[i][k_idx][z_0-1] + (*a)[i][k_idx][z_0]);
  return a_z0;
}

mat_cmplx average_array_z0(const tens_cmplx &a, const unsigned z_0)
{
  using namespace parameters::numerical;

  mat_cmplx a_z0(a.size(), vec_cmplx(a[0].size(), 0.0));

  for (unsigned i = 0; i < a.size(); ++i)
    for (unsigned k_idx = 0; k_idx < a[i].size(); ++k_idx)
      a_z0[i][k_idx] = 0.5 * (a[i][k_idx][z_0-1] + a[i][k_idx][z_0]);
  return a_z0;
}


// The accuracy is measured by the largest relative deviation of the a_i's with respect to the previous iteration
double update_accuracy_z_dep(const dressing* a, const dressing* a_old)
{
  double current_acc = 0.;

  using namespace parameters::numerical;
  for (unsigned i = 0; i < n_structs; ++i)
    for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
    {
      double dif = 0.;
      for (unsigned z_idx = 0; z_idx < (*a)[i][k_idx].size(); ++z_idx)
      {
        const auto sum = (*a)[i][k_idx][z_idx] + (*a_old)[i][k_idx][z_idx];
        dif += std::abs((*a)[i][k_idx][z_idx] - (*a_old)[i][k_idx][z_idx]) / double((*a)[i][k_idx].size());
        if(!isEqual(std::abs(sum), 0.))
          current_acc = std::max(current_acc, dif / std::abs(sum) );        
      }
    }
  return current_acc;
}


void a_iteration_step(
    dressing* a, const dressing* b, BSE_kernel_L* K_prime_L, BSE_kernel_T* K_prime_T, // BSE dressings, b and K^\prime_ij
    const vec_double &z_grid, const vec_double &log_k_sq_grid, // angular and momentum variable grids
    const Integrator2d& qint2d, const DiscrIntegrator2d& qy_int2d,// 2D integration
    const double z2 // quark renormalization constant
    )
{
    using namespace parameters::numerical;

    const vec_cmplx temp0(z_steps, 0.0);
    const mat_cmplx temp1(k_steps, temp0);

    // calculate transverse dressing functions

    for (unsigned i = 0; i < n_structs_T; ++i) // loop over tensor structures
    #pragma omp parallel for collapse(2)
        for (unsigned i_k = 0; i_k < k_steps; ++i_k) // loop over k grid
            for (unsigned i_z = 0; i_z < z_steps; ++i_z) // loop over z grid
            {
              // Initialize a_i with inhomogeneous term
              (*a)[i][i_k][i_z] = z2 * a0(i);

              // initialize 2d function which will be integrated later
              const vec_cmplx temp1(z_steps,0);
              mat_cmplx b_K_prime(k_steps,temp1);

              for (unsigned j =0; j < n_structs_T; ++j)
              {
                // If kernel matrix element vanishes, skip integration
                if (K::isZeroIndex(i, j))
                    continue;
                
                for(unsigned k_prime_idx = 0; k_prime_idx < k_steps; ++k_prime_idx)
                  for(unsigned z_prime_idx = 0; z_prime_idx < z_steps; ++z_prime_idx)
                  {
                    const double& z_prime = z_grid[z_prime_idx];
                    const double log_k_prime_sq = log_k_sq_grid[k_prime_idx];

                    b_K_prime[k_prime_idx][z_prime_idx] += jakobian_log_k_prime(z_prime, log_k_prime_sq) *
                      (*b)[j][k_prime_idx][z_prime_idx] * (*K_prime_T)[i][i_k][i_z][j][k_prime_idx][z_prime_idx];               
                  }
              }

              // Perform integral
              const std::complex<double> integral = qy_int2d(b_K_prime,log_k_sq_grid[0], log_k_sq_grid[k_steps - 1],z_grid[0], z_grid[z_steps - 1]);

              // Add integral to initialized value of a_i
              (*a)[i][i_k][i_z] += integral;
            }


    // calculate longitudinal dressing functions

    for (unsigned i = n_structs_T; i < n_structs; ++i) // loop over tensor structures
    #pragma omp parallel for collapse(2)
        for (unsigned i_k = 0; i_k < k_steps; ++i_k) // loop over k grid
            for (unsigned i_z = 0; i_z < z_steps; ++i_z) // loop over z grid
            {
              // Initialize a_i with inhomogeneous term
              (*a)[i][i_k][i_z] = z2 * a0(i);

              // initialize 2d function which will be integrated later
              const vec_cmplx temp1(z_steps,0);
              mat_cmplx b_K_prime(k_steps,temp1);

              for (unsigned j = n_structs_T; j < n_structs; ++j)
              {
                // If kernel matrix element vanishes, skip integration
                if (K::isZeroIndex(i, j))
                    continue;
                
                for(unsigned k_prime_idx = 0; k_prime_idx < k_steps; ++k_prime_idx)
                  for(unsigned z_prime_idx = 0; z_prime_idx < z_steps; ++z_prime_idx)
                  {
                    const double& z_prime = z_grid[z_prime_idx];
                    const double log_k_prime_sq = log_k_sq_grid[k_prime_idx];

                    b_K_prime[k_prime_idx][z_prime_idx] += jakobian_log_k_prime(z_prime, log_k_prime_sq) *
                      (*b)[j][k_prime_idx][z_prime_idx] * (*K_prime_L)[i-8][i_k][i_z][j-8][k_prime_idx][z_prime_idx];               
                  }
              }

              // Perform integral
              const std::complex<double> integral = qy_int2d(b_K_prime,log_k_sq_grid[0], log_k_sq_grid[k_steps - 1],z_grid[0], z_grid[z_steps - 1]);

              // Add integral to initialized value of a_i
              (*a)[i][i_k][i_z] += integral;
            }
}


void b_iteration_step(const dressing* a, const double &q_sq, const vec_double &z_grid, const vec_double &log_k_sq_grid, dressing* b, propagator_kernel* kernel_G)
{
  using namespace parameters::numerical;
  #pragma omp parallel for collapse(2) // parallelize the outermost two loops
  for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
    for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx) 
      for (unsigned i_global = 0; i_global < n_structs; ++i_global)
      {
        // Initialize the b's to 0
        (*b)[i_global][k_idx][z_idx] = 0.0;

        unsigned i_block = floor(i_global / 4.);
      
        for (unsigned j_local = 0; j_local < 4; ++j_local)
        {
          // since the global iterators i,j run from 0 to 11, we need to map them to the corresponding local indices of the propagator matrix
          unsigned i_local = i_global - 4 * i_block;
          unsigned j_global = j_local + 4 * i_block;

          // perform multiplication of propagator matrix G and BSE dressing vector a
          (*b)[i_global][k_idx][z_idx] += (*kernel_G)[i_block][i_local][j_local][k_idx][z_idx] * (*a)[j_global][k_idx][z_idx];
        }
      }
}


void calculate_fg(dressing* a, dressing* fg, const double& q_sq, const vec_double &log_k_sq_grid, const vec_double &z_grid)
{
  using namespace parameters::numerical;
  using namespace basistransform;

  vec_cmplx a_copy(n_structs);

  const double Q = std::sqrt(q_sq);

  for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
  {
    double k_sq = std::exp(log_k_sq_grid[k_idx]);
    double k = std::sqrt(k_sq);

    for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx)
    {
      for (unsigned i = 8; i < n_structs; ++i)
      {  // evaluate a at given k_sq and q_sq
        a_copy[i] = (*a)[i][k_idx][z_idx];
      }

      const double z = z_grid[z_idx];
      const double s = std::sqrt(1. - powr<2>(z));

      (*fg)[0][k_idx][z_idx] = f1(Q, s, z, k, a_copy);
      (*fg)[1][k_idx][z_idx] = f2(Q, s, z, k, a_copy);
      (*fg)[2][k_idx][z_idx] = f3(Q, s, z, k, a_copy);
      (*fg)[3][k_idx][z_idx] = f4(Q, s, z, k, a_copy);
      (*fg)[4][k_idx][z_idx] = f5(Q, s, z, k, a_copy);
      (*fg)[5][k_idx][z_idx] = f6(Q, s, z, k, a_copy);
      (*fg)[6][k_idx][z_idx] = f7(Q, s, z, k, a_copy);
      (*fg)[7][k_idx][z_idx] = f8(Q, s, z, k, a_copy);
      (*fg)[8][k_idx][z_idx] = g1(Q, s, z, k, a_copy);
      (*fg)[9][k_idx][z_idx] = g2(Q, s, z, k, a_copy);
      (*fg)[10][k_idx][z_idx] = g3(Q, s, z, k, a_copy);
      (*fg)[11][k_idx][z_idx] = g4(Q, s, z, k, a_copy);
    }
  }
}


template<typename Quark>
void solve_BSE(const vec_double &q_grid, const vec_double &z_grid, const vec_double &log_k_sq_grid, const vec_double &y_grid, const bool use_PauliVillars, const bool debug)
{
  using namespace std::chrono;
  const auto start_time = steady_clock::now();

  using namespace parameters::numerical;
  const unsigned z_0 = z_grid.size() / 2;


  // initialze some complex vectors and matrices for later use
  const vec_cmplx temp0(z_steps, 0.0);
  const mat_cmplx temp1(k_steps, temp0);

   
  // Do some Legendre Magic
  DiscrIntegrator1d yint1d;
  DiscrIntegrator2d qy_int2d;
  Integrator1d qint1d;
  Integrator2d qint2d;


  // initialize quark propagator. if the DSE is used, this automatically solves the DSE.
  std::cout << "\n\n_______________    Quark DSE    _______________\n\n\n";
  const Quark quark;  
  const auto DSE_end_time = steady_clock::now();
  std::cout << ", t = " << duration_cast<milliseconds>(DSE_end_time - start_time).count()/1000.<< "s)\n";


  // prepare output files
  emptyIdxFile<12>("a_file", "#q_sq i k_sq z Re(fg) Im(fg)");
  emptyIdxFile<12>("fg_file", "#q_sq i k_sq z Re(fg) Im(fg)");
  emptyIdxFile<12>("fg_z0_file", "#q_sq i k_sq Re(fg) Im(fg)");
  emptyIdxFile<3>("w_file", "#q_sq i k_sq z Re(w) Im(w)");
  emptyIdxFile<3>("w_z0_file", "#q_sq i k_sq Re(w) Im(w)");

  std::cout << "\nCalculating the WTIs..." << std::flush;

  // check if quark DSE or quark model is used
  string quarktype = type_name<Quark>();
  
  vec_double quark_M(quark.quark_grid.size());

  if (quarktype.compare("quark_DSE") == 0) {

    // construct mass function from quark dressings
    for (unsigned i = 0; i < quark_M.size(); ++i)
      quark_M[i] = quark.quark_b[i] / quark.quark_a[i];   
  }

  // construct interpolating functions for quark dressing and mass function
  tk::spline ip_A(quark.quark_grid, quark.quark_a);
  tk::spline ip_M(quark.quark_grid, quark_M);

  // save WTI
  for (unsigned q_iter = 0; q_iter < q_steps; q_iter++)
  {
    const double q_sq = q_grid[q_iter];
    tens_cmplx w(3, temp1);

    for (unsigned k_idx = 0; k_idx < parameters::numerical::k_steps; ++k_idx)
      for (unsigned z_idx = 0; z_idx < parameters::numerical::z_steps; ++z_idx)
      {
        const double k_sq = std::exp(log_k_sq_grid[k_idx]);
        const double& z = z_grid[z_idx];

        double kplus2 = k_sq + q_sq/4. + std::sqrt(q_sq * k_sq) * z;
        double kminus2 = k_sq + q_sq/4. - std::sqrt(q_sq * k_sq) * z;
        double quark_A_p, quark_A_m, quark_M_p, quark_M_m;

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        // OUTSOURCE CONSTRUCTION OF INTERPOLATING FUNCTION TO A CLASS 

        if (quarktype.compare("quark_DSE") == 0) {          
          quark_A_p = ip_A(std::log(kplus2));
          quark_A_m = ip_A(std::log(kminus2));
          quark_M_p = ip_M(std::log(kplus2));
          quark_M_m = ip_M(std::log(kminus2));
        } else {
          quark_A_p = quark.A(std::log(kplus2));
          quark_A_m = quark.A(std::log(kminus2));
          quark_M_p = quark.M(std::log(kplus2));
          quark_M_m = quark.M(std::log(kminus2));
        };

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        w[0][k_idx][z_idx] = Sigma_A(kminus2,kplus2,quark_M_p,quark_M_m,quark_A_p,quark_A_m);
        w[1][k_idx][z_idx] = Delta_A(kminus2,kplus2,quark_M_p,quark_M_m,quark_A_p,quark_A_m);
        w[2][k_idx][z_idx] = Delta_B(kminus2,kplus2,quark_M_p,quark_M_m,quark_A_p,quark_A_m);
      }

    saveToFile_withGrids<3>(w, "w_file", q_sq, log_k_sq_grid, z_grid);
    const auto w_z0 = average_array_z0(w, z_0);
    saveToFile_withGrids<3>(w_z0, "w_z0_file", q_sq, log_k_sq_grid);
  }

  std::cout << " done\n";

std::cout << "\n\n_______________    Quark-photon-vertex BSE    _______________\n\n\n";

  // Calculate BSE kernel
    std::cout << "\nCalculating BSE kernel..." << std::flush;

    const auto BSE_kernel_start_time = steady_clock::now();

    // allocate
    BSE_kernel_L *kernel_K_L = new BSE_kernel_L; // allocate
    BSE_kernel_T *kernel_K_T = new BSE_kernel_T;

    calculate_BSE_kernel(y_grid, yint1d, z_grid, log_k_sq_grid, kernel_K_L, kernel_K_T, quark.z2(), use_PauliVillars);

    // end BSE kernel timer
    const auto BSE_kernel_end_time = steady_clock::now();
    std::cout << " done (t = " << duration_cast<milliseconds>(BSE_kernel_end_time - BSE_kernel_start_time).count()/1000.<< "s)\n";

  // loop over q
  for (unsigned q_iter = 0; q_iter < q_grid.size(); q_iter++)
  {
    // start iteration timer
    const auto iter_start_time = steady_clock::now();

    const double &q_sq = q_grid[q_iter];

    std::cout <<  "\n\nCalculation for q^2 = " << q_sq << "\n";

    // initialize a, b, fg
    dressing *a = new dressing;
    dressing *b = new dressing;
    dressing *fg = new dressing;

    // Calculate propagator kernel
    std::cout << " - Calculating propagator kernel...";

    propagator_kernel *kernel_G = new propagator_kernel; // allocate

    // calculate_propagator_kernel<Quark>(q_sq, z_grid, log_k_sq_grid, quark, kernel_G);
    calculate_propagator_kernel<Quark>(q_sq, z_grid, log_k_sq_grid, quark, kernel_G, ip_A, ip_M);

    // end propagator kernel timer
    std::cout << " done\n";

    // Initialize a with bare vertex
    a_initialize(a, quark.z2());

    std::cout << " - Solving equation..." << std::flush;
    double current_acc = 1.0;
    unsigned current_step = 0;
    while (max_steps > current_step++ && current_acc > target_acc)
    {
      debug_out("\n    Started a step...\n", debug);
      
      // copy for checking the convergence
      dressing* a_old = new dressing;
      *a_old = *a;
      
      debug_out("    Calculating b_i...", debug);
      b_iteration_step(a, q_sq, z_grid, log_k_sq_grid, b, kernel_G);
      debug_out(" done\n", debug);

      debug_out("    Calculating a_i...", debug);
      a_iteration_step(a, b, kernel_K_L, kernel_K_T, z_grid, log_k_sq_grid, qint2d, qy_int2d, quark.z2());
      debug_out(" done\n", debug);

      // check the convergence
      current_acc = update_accuracy_z_dep(a, a_old);
      debug_out("    current_step = " + std::to_string(current_step) + "\n    current_acc = " + std::to_string(current_acc) + "\n", debug);

      delete a_old;
    }

    // stop the iteration timer
    const auto iter_end_time = steady_clock::now();

    if (current_acc < target_acc)
      std::cout << " done (" << current_step << " iterations, t = " << duration_cast<milliseconds>(iter_end_time - iter_start_time).count()/1000.<< "s)\n";
    else
      std::cout << "  ! Did not converge !\n";

    saveToFile_withGrids<n_structs>(a, "a_file", q_sq, log_k_sq_grid, z_grid);

    std::cout << " - Saving results..." << std::flush;
    // transform to g,f (almost in place!)
    calculate_fg(a, fg, q_sq, log_k_sq_grid, z_grid);

    // save to the prepared file
    saveToFile_withGrids<n_structs>(fg, "fg_file", q_sq, log_k_sq_grid, z_grid);

    // calculate angular average
    const auto fg_z0 = average_array_z0(fg, z_0);
    saveToFile_withGrids<n_structs>(fg_z0, "fg_z0_file", q_sq, log_k_sq_grid);
    std::cout << "  done\n";

    std::cout << "Calculation finished after " << duration_cast<milliseconds>(iter_end_time - iter_start_time).count()/1000.<< "s\n";
    
    delete a;
    delete b;
    delete fg;
    delete kernel_G;
  }

  delete kernel_K_L;
  delete kernel_K_T;

  auto end_time = steady_clock::now();
  std::cout << "\nProgram finished after " << duration_cast<milliseconds>(end_time - start_time).count()/1000.<< "s\n";
}


template<typename Quark>
void solve_DSE(const vec_double &q_grid, const vec_double &z_grid, const vec_double &log_k_sq_grid, const vec_double &y_grid, const bool use_PauliVillars, const bool debug)
{
  using namespace parameters::numerical;
  const unsigned z_0 = 1;

  using namespace std::chrono;
  const auto start_time = steady_clock::now();

  const vec_cmplx temp0(parameters::numerical::z_steps, 0.0);
  const mat_cmplx temp1(parameters::numerical::k_steps, temp0);

  const Quark quark;

  // prepare output files
  emptyIdxFile<3>("w_z0_file", "#q_sq i k_sq Re(w) Im(w)");
  emptyFile("quark_M_Z", "#k_sq M Z");
  emptyFile("test", "abc");

  vec_double quark_M(parameters::numerical::k_steps,0.0);
  vec_double quark_Z(parameters::numerical::k_steps,0.0);

    std::cout << "\nCalculating the WTIs...\n";

  // check WTI
  for (unsigned q_iter = 0; q_iter < q_steps; q_iter++)
  {
    double q_sq = q_grid[q_iter];
    tens_cmplx w(3, temp1);
    for (unsigned k_idx = 0; k_idx < parameters::numerical::k_steps; ++k_idx)
    {
      const double k_sq = std::exp(log_k_sq_grid[k_idx]);
      for (unsigned z_idx = 0; z_idx < parameters::numerical::z_steps; ++z_idx)
      {
        const double& z = z_grid[z_idx];

        double kplus2 = k_sq + q_sq/4. + std::sqrt(q_sq * k_sq) * z;
        double kminus2 = k_sq + q_sq/4. - std::sqrt(q_sq * k_sq) * z;

        w[0][k_idx][z_idx] = Sigma_A<Quark>(kminus2,kplus2,quark);
        w[1][k_idx][z_idx] = Delta_A<Quark>(kminus2,kplus2,quark);
        w[2][k_idx][z_idx] = Delta_B<Quark>(kminus2,kplus2,quark);
      }
      quark_M[k_idx] = quark.M(k_sq);
      quark_Z[k_idx] = 1. / quark.A(k_sq);
    }

    const auto w_z0 = average_array_z0(w, z_0);
    std::cout << "w_0(10^-4) = " << w_z0[0][0] << std::endl;
    std::cout << "w_1(10^-4) = " << w_z0[1][0] << std::endl;
    std::cout << "w_2(10^-4) = " << w_z0[2][0] << std::endl;
    saveToFile_withGrids<3>(w_z0, "w_z0_file", q_sq, log_k_sq_grid);
    saveToFile_M_Z(quark_M, quark_Z, log_k_sq_grid, "quark_M_Z");
  }
} 