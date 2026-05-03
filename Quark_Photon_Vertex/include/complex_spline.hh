#pragma once

#include <complex>
#include <vector>

#include "spline.h"
#include "types.hh"

// Cubic spline of a complex-valued grid function in one real coordinate.
// Two natural cubic splines (re/im), sharing the x grid. Used to lift the
// log-k' interpolation in BSE/HVP integrands from linear to cubic — linear
// was the dominant residual at the IR k corner where WTI2 max-error sits.
struct ComplexSpline {
  tk::spline re_, im_;

  void set_points(const std::vector<double>& x, const vec_cmplx& y) {
    std::vector<double> y_re(y.size()), y_im(y.size());
    for (std::size_t i = 0; i < y.size(); ++i) {
      y_re[i] = y[i].real();
      y_im[i] = y[i].imag();
    }
    re_.set_points(x, y_re);
    im_.set_points(x, y_im);
  }

  std::complex<double> operator()(double x) const {
    return std::complex<double>(re_(x), im_(x));
  }
};
