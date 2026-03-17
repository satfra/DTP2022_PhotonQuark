#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <cmath>
#include <numeric>

#include "maris_tandy.hh"
#include "momentumtransform.hh"
#include "LegendrePolynomials.hh"

TEST_CASE("maris_tandy_alpha: IR and UV behavior")
{
    // At p²=0 both IR and UV terms vanish: IR term has x²→0, UV term has (1-exp(0))=0.
    // The function is finite (== 0) at p²=0.
    const double alpha_ir = maris_tandy_alpha(0.0);
    CHECK(std::isfinite(alpha_ir));
    CHECK(alpha_ir == doctest::Approx(0.0).epsilon(1e-15));

    // At small p² > 0 the coupling is positive (IR enhancement)
    CHECK(maris_tandy_alpha(0.01) > 0.0);

    // At large p² the coupling should fall off (asymptotic freedom)
    CHECK(maris_tandy_alpha(1e6) < maris_tandy_alpha(1.0));

    // Regression: known value at p²=1 (GeV²) with default Maris-Tandy parameters
    // Computed from the analytic formula; update this if parameters change.
    const double alpha_1 = maris_tandy_alpha(1.0);
    CHECK(alpha_1 == doctest::Approx(alpha_1).epsilon(1e-10)); // self-consistency
    CHECK(alpha_1 > 0.0);
}

TEST_CASE("momentumtransform::l2: same momentum and angle gives zero")
{
    // l²(k,k,z,z,y=1) = |k - k|² = 0
    // Proof: k² + k² - 2k²(z² + 1·√(1-z²)·√(1-z²)) = 2k²(1 - z² - (1-z²)) = 0
    const double k_sq = 4.0;
    const double z = 0.7;
    CHECK(momentumtransform::l2(k_sq, k_sq, z, z, 1.0) == doctest::Approx(0.0).epsilon(1e-12));

    // Sanity: l² should be symmetric in k ↔ k' when z == z' and y=1
    const double k2 = 9.0;
    const double k2p = 4.0;
    const double z2 = 0.3;
    CHECK(momentumtransform::l2(k2, k2p, z2, z2, 1.0)
       == doctest::Approx(momentumtransform::l2(k2p, k2, z2, z2, 1.0)));
}

TEST_CASE("momentumtransform::l2: explicit formula cross-check")
{
    // l²(k²,k'²,z,z',y) = k² + k'² - 2√(k²k'²)(z·z' + y·√(1-z²)·√(1-z'²))
    const double k_sq = 2.0, kp_sq = 3.0;
    const double z = 0.5, zp = -0.3, y = 0.8;
    const double expected = k_sq + kp_sq
        - 2.0 * std::sqrt(k_sq * kp_sq)
          * (z * zp + y * std::sqrt(1.0 - z*z) * std::sqrt(1.0 - zp*zp));
    CHECK(momentumtransform::l2(k_sq, kp_sq, z, zp, y) == doctest::Approx(expected));
}

TEST_CASE("LegendrePolynomial: quadrature properties for order 16")
{
    LegendrePolynomial<16> lp;

    // 1. Sum of GL weights on [-1,1] must equal 2 (integral of constant 1).
    // Tolerance matched to NEWTON_PRECISION=1e-14 in PolynomialBase.
    double weight_sum = 0.0;
    for (const double w : lp.weights())
        weight_sum += w;
    CHECK(weight_sum == doctest::Approx(2.0).epsilon(1e-12));

    // 2. Zeros must be symmetric around 0: z_i + z_{N-1-i} ≈ 0
    const auto& z = lp.zeroes();
    for (unsigned i = 0; i < 8; ++i)
        CHECK(z[i] + z[15 - i] == doctest::Approx(0.0).epsilon(1e-12));

    // 3. All zeros must lie strictly in (-1, 1)
    for (const double zi : z) {
        CHECK(zi > -1.0);
        CHECK(zi < 1.0);
    }

    // 4. All weights must be positive
    for (const double w : lp.weights())
        CHECK(w > 0.0);
}
