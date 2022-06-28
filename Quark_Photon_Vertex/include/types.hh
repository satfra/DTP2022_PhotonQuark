#pragma once

#include <complex>
#include <vector>
#include "parameters.hh"
#include <array>

using namespace::parameters::numerical;

using vec_cmplx = std::vector<std::complex<double>>;
using mat_cmplx = std::vector<vec_cmplx>;
using tens_cmplx = std::vector<mat_cmplx>;
using tens2_cmplx = std::vector<tens_cmplx>;
using jtens2_cmplx = std::vector<tens2_cmplx>;
using ijtens2_cmplx = std::vector<jtens2_cmplx>;

using vec_double = std::vector<double>;
using mat_double = std::vector<vec_double>;
using tens_double = std::vector<mat_double>;
using tens2_double = std::vector<tens_double>;
using jtens2_double = std::vector<tens2_double>;
using ijtens2_double = std::vector<jtens2_double>;

using dressing = std::array<std::array <std::array <std::complex<double>, z_steps>, k_steps>, n_structs>;
using WTI = std::array<std::array <std::array <std::complex<double>, z_steps>, k_steps>, 3>;
using BSE_kernel_T = std::array<std::array <std::array <std::array <std::array <std::array <double, z_steps>, k_steps>, n_structs_T>, z_steps>, k_steps>, n_structs_T>;
using BSE_kernel_L = std::array<std::array <std::array <std::array <std::array <std::array <double, z_steps>, k_steps>, n_structs_L>, z_steps>, k_steps>, n_structs_L>;
using propagator_kernel = std::array< std::array< std::array <std::array <std::array <std::complex<double>, z_steps>, k_steps>, n_structs_L>, n_structs_L>, 3>;
