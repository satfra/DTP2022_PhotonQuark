#pragma once

#include <cmath>
#include <Utils.hh>
#include <PolynomialBase.hh>

// Gauss-Chebyshev quadrature of the SECOND kind (polynomials U_n).
// Approximates  ∫_{-1}^{1} f(z) · √(1 − z²) dz  ≈  Σ_k w_k · f(z_k)
// with closed-form nodes/weights:
//   z_k = cos(k π / (n+1))
//   w_k = π/(n+1) · sin²(k π / (n+1))
// Exact for polynomial f of degree ≤ 2n − 1.
//
// IMPORTANT: When used with qIntegral / qIntegral2d, the integration
// bounds MUST be (-1, 1). Linear remapping to [a, b] is correct for
// Legendre (constant weight) but breaks the meaning of the implicit
// √(1 − z²) weight here.
template<unsigned o, typename _RF = double>
class ChebyshevPolynomial2 : public PolynomialBase<ChebyshevPolynomial2<o>, _RF>
{
  using PB = PolynomialBase<ChebyshevPolynomial2<o>, _RF>;

  public:
    static constexpr unsigned order = o;
    using RF = _RF;

    ChebyshevPolynomial2()
    {
      calcZeros();
      calcWeights();
    }

  protected:
    // U_n via the same recurrence as T_n but with U_1 = 2x (vs T_1 = x).
    virtual RF P(const RF& x) const override
    {
      std::vector<RF> res(order+1);
      res[0] = 1.;
      res[1] = 2.*x;
      for(unsigned n = 2; n <= order; ++n)
        res[n] = 2.*x * res[n-1] - res[n-2];
      return res[order];
    }

    virtual RF dP(const RF& x) const override
    {
      std::vector<RF> der(order+1);
      std::vector<RF> res(order+1);
      res[0] = 1.;
      der[0] = 0.;
      res[1] = 2.*x;
      der[1] = 2.;
      for(unsigned n = 2; n <= order; ++n)
      {
        res[n] = 2.*x * res[n-1] - res[n-2];
        der[n] = 2.* res[n-1] + 2.*x * der[n-1] - der[n-2];
      }
      return der[order];
    }

    // Never called: calcZeros() below uses the closed form.
    virtual RF initialGuess(unsigned j) const override
    {
      return std::cos(M_PI * RF(j+1) / (order+1));
    }

    // Closed-form zeros, ascending. Index i ↔ k = order − i so that
    // z_0 = cos(order·π/(order+1))  (most negative, near −1)
    // z_{order-1} = cos(π/(order+1)) (most positive, near +1)
    // — matches the ascending ordering produced by Newton's method on
    // LegendrePolynomial.
    void calcZeros() override
    {
      PB::_z.resize(order);
      for(unsigned i = 0; i < order; ++i)
      {
        const unsigned k = order - i;
        PB::_z[i] = std::cos(M_PI * RF(k) / RF(order + 1));
      }
    }

    void calcWeights() override
    {
      PB::_w.resize(order);
      for(unsigned i = 0; i < order; ++i)
      {
        const unsigned k = order - i;
        const RF s = std::sin(M_PI * RF(k) / RF(order + 1));
        PB::_w[i] = M_PI / RF(order + 1) * s * s;
      }
    }
};
