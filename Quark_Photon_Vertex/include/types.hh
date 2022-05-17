#pragma once

#include <complex>
#include <vector>

using vec_cmplx = std::vector<std::complex<double>>;
using mat_cmplx = std::vector<vec_cmplx>;
using tens_cmplx = std::vector<mat_cmplx>;

using vec_double = std::vector<double>;
using mat_double = std::vector<vec_double>;
using tens_double = std::vector<mat_double>;
using tens2_double = std::vector<tens_double>;
using jtens2_double = std::vector<tens2_double>;
using ijtens2_double = std::vector<jtens2_double>;
