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
#include <chrono>
//#include "quark_dressings.hh"

using namespace std;
using namespace parameters::numerical;

using Integrator1d = qIntegral<LegendrePolynomial<parameters::numerical::y_steps>>;
using DiscrIntegrator1d = yIntegral<LegendrePolynomial<parameters::numerical::y_steps>>;
using Integrator2d = qIntegral2d<LegendrePolynomial<parameters::numerical::k_steps>, LegendrePolynomial<parameters::numerical::z_steps>>;
using DiscrIntegrator2d = gaulegIntegral2d<LegendrePolynomial<parameters::numerical::k_steps>, LegendrePolynomial<parameters::numerical::z_steps>>;



template<typename Quark>
void calculate_propagator_kernel(const double &q_sq, const vec_double &z_grid, const vec_double &log_k_sq_grid, const Quark& quark, propagator_kernel* kernel_G, const tk::spline& ip_A, const tk::spline& ip_M)
{
  using namespace parameters::numerical;

  #pragma omp parallel for collapse(2) // parallelize the outermost two loops
  for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
    for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx) 
    {
      const double k_sq = std::exp(log_k_sq_grid[k_idx]);
      const double& z = z_grid[z_idx];

      const double k_p_sq = std::log(k_sq + 0.25 * q_sq + std::sqrt(k_sq * q_sq) * z);
      const double k_m_sq = std::log(k_sq + 0.25 * q_sq - std::sqrt(k_sq * q_sq) * z);

      double quark_A_p, quark_A_m, quark_M_m, quark_M_p;

      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      // OUTSOURCE CONSTRUCTION OF INTERPOLATING FUNCTION TO A CLASS 

      string quarktype = type_name<Quark>();

      if (quarktype.compare("quark_DSE") == 0) {   
        quark_A_p = ip_A(k_p_sq);
        quark_A_m = ip_A(k_m_sq);
        quark_M_p = ip_M(k_p_sq);
        quark_M_m = ip_M(k_m_sq);
      } else {
        quark_A_p = quark.A(k_p_sq);
        quark_A_m = quark.A(k_m_sq);
        quark_M_p = quark.M(k_p_sq);
        quark_M_m = quark.M(k_m_sq);
      };

      // quark_dressings<Quark> _quark_dressings(quark);
      // quark_A_p = _quark_dressings.A(k_p_sq);
      // quark_A_m = _quark_dressings.A(k_m_sq);
      // quark_M_p = _quark_dressings.M(k_p_sq);
      // quark_M_m = _quark_dressings.M(k_m_sq);

      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      // initialization of propagator kernel. this already calculates all the entries for a given kinematic configuration 
      const G G_(k_sq, z, q_sq, quark_M_p, quark_M_m, quark_A_p, quark_A_m);

      
      for (unsigned i_global = 0; i_global < n_structs; ++i_global)
      {
        unsigned i_block = (i_global == 0) ? 0 : floor(i_global / 4.);
       
        for (unsigned j_local = 0; j_local < 4; ++j_local)
        {
          // since the global iterators i,j run from 0 to 11, we need to map them to the corresponding local indices of the propagator matrix
          unsigned i_local = i_global - 4 * i_block;
          unsigned j_global = j_local + 4 * i_block;
          
          // evaluate kernel and store
          (*kernel_G)[i_block][i_local][j_local][k_idx][z_idx] = G_.get(i_global, j_global);
        }
      }
    }  
}
