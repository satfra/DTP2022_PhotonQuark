#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <cmath>

#include "Utils.hh"
#include "LegendrePolynomials.hh"
#include "ChebyshevPolynomial2.hh"
#include "QuadratureIntegral.hh"

TEST_CASE("powr: compile-time integer power")
{
    CHECK(powr<0>(5.0) == 1.0);
    CHECK(powr<0>(0.0) == 1.0);
    CHECK(powr<1>(3.0) == doctest::Approx(3.0));
    CHECK(powr<2>(3.0) == doctest::Approx(9.0));
    CHECK(powr<3>(2.0) == doctest::Approx(8.0));
    CHECK(powr<-1>(2.0) == doctest::Approx(0.5));
    CHECK(powr<-2>(2.0) == doctest::Approx(0.25));
}

TEST_CASE("isEqual: floating-point equality")
{
    CHECK(isEqual(1.0, 1.0));
    CHECK(isEqual(0.0, 0.0));
    CHECK(!isEqual(1.0, 2.0));

    // values within machine epsilon are equal
    const double eps = std::numeric_limits<double>::epsilon();
    CHECK(isEqual(1.0, 1.0 + eps));

    // values differing by 1 ULP at large scale are NOT equal (1/1e10 >> eps)
    CHECK(!isEqual(1e10, 1e10 + 1.0));

    // custom tolerance
    CHECK(isEqual(1.0, 1.05, 0.1));
    CHECK(!isEqual(1.0, 1.2, 0.1));
}

TEST_CASE("linearMapTo: scalar version")
{
    // identity map
    CHECK(linearMapTo(0.5, 0.0, 1.0, 0.0, 1.0) == doctest::Approx(0.5));

    // reverse map
    CHECK(linearMapTo(0.5, 0.0, 1.0, 1.0, 0.0) == doctest::Approx(0.5));

    // boundary values
    CHECK(linearMapTo(0.0, 0.0, 1.0, -1.0, 1.0) == doctest::Approx(-1.0));
    CHECK(linearMapTo(1.0, 0.0, 1.0, -1.0, 1.0) == doctest::Approx(1.0));

    // midpoint maps to midpoint
    CHECK(linearMapTo(0.5, 0.0, 1.0, -1.0, 1.0) == doctest::Approx(0.0));

    // arbitrary range
    CHECK(linearMapTo(5.0, 0.0, 10.0, 100.0, 200.0) == doctest::Approx(150.0));
}

TEST_CASE("linearMapTo: vector version")
{
    std::vector<double> v = {0.0, 0.5, 1.0};
    auto mapped = linearMapTo(v, 0.0, 1.0, -1.0, 1.0);
    REQUIRE(mapped.size() == 3);
    CHECK(mapped[0] == doctest::Approx(-1.0));
    CHECK(mapped[1] == doctest::Approx(0.0));
    CHECK(mapped[2] == doctest::Approx(1.0));

    // original vector is unchanged
    CHECK(v[0] == 0.0);
}

TEST_CASE("locate: binary search")
{
    std::vector<double> xx = {0.0, 1.0, 2.0, 3.0, 4.0};

    // middle value
    CHECK(locate(xx, 1.5) == 1);

    // left boundary: locate returns 0 for xx[0]
    CHECK(locate(xx, 0.0) == 0);

    // right boundary: locate returns n-2 for xx[n-1]
    CHECK(locate(xx, 4.0) == 3);

    // descending vector
    std::vector<double> desc = {4.0, 3.0, 2.0, 1.0, 0.0};
    const int j = locate(desc, 2.5);
    CHECK(j >= 0);
    CHECK(j < (int)desc.size() - 1);
    CHECK((desc[j] >= 2.5 || desc[j+1] >= 2.5)); // x is between xx[j] and xx[j+1]
}

TEST_CASE("ChebyshevPolynomial2: closed-form sanity")
{
    ChebyshevPolynomial2<16> cp;
    const auto& z = cp.zeroes();
    const auto& w = cp.weights();

    REQUIRE(z.size() == 16);
    REQUIRE(w.size() == 16);

    // All zeros strictly inside (-1, 1)
    for (const double zi : z) {
        CHECK(zi > -1.0);
        CHECK(zi < 1.0);
    }

    // Ascending order (matches LegendrePolynomial convention)
    for (unsigned i = 1; i < z.size(); ++i)
        CHECK(z[i] > z[i-1]);

    // Symmetric about 0
    for (unsigned i = 0; i < 8; ++i)
        CHECK(z[i] + z[15 - i] == doctest::Approx(0.0).epsilon(1e-14));

    // All weights positive
    for (const double wi : w)
        CHECK(wi > 0.0);
}

TEST_CASE("ChebyshevPolynomial2: weight sum equals pi/2 = ∫_{-1}^{1} √(1-z²) dz")
{
    ChebyshevPolynomial2<16> cp;
    double s = 0.0;
    for (const double wi : cp.weights()) s += wi;
    CHECK(s == doctest::Approx(M_PI_2).epsilon(1e-14));
}

TEST_CASE("ChebyshevPolynomial2: polynomial moments against √(1-z²)")
{
    // ∫_{-1}^{1} z^(2k) √(1-z²) dz = π/2 · (2k)!/(4^k · k! · (k+1)!)
    // Concretely:  k=0 → π/2 ; k=1 → π/8 ; k=2 → π/16 ; k=3 → 5π/128
    qIntegral<ChebyshevPolynomial2<16>> qint;

    auto m0 = [](const double& z) { (void)z; return 1.0; };
    auto m2 = [](const double& z) { return z*z; };
    auto m4 = [](const double& z) { return z*z*z*z; };
    auto m6 = [](const double& z) { return z*z*z*z*z*z; };
    auto m1 = [](const double& z) { return z; };
    auto m3 = [](const double& z) { return z*z*z; };

    CHECK(qint(m0, -1.0, 1.0) == doctest::Approx(M_PI / 2.0).epsilon(1e-14));
    CHECK(qint(m2, -1.0, 1.0) == doctest::Approx(M_PI / 8.0).epsilon(1e-14));
    CHECK(qint(m4, -1.0, 1.0) == doctest::Approx(M_PI / 16.0).epsilon(1e-14));
    CHECK(qint(m6, -1.0, 1.0) == doctest::Approx(5.0 * M_PI / 128.0).epsilon(1e-14));

    // odd moments vanish
    CHECK(std::abs(qint(m1, -1.0, 1.0)) < 1e-14);
    CHECK(std::abs(qint(m3, -1.0, 1.0)) < 1e-14);
}

TEST_CASE("ChebyshevPolynomial2: cross-check vs Legendre on f(z)·√(1-z²)")
{
    // ∫_{-1}^{1} exp(z) √(1-z²) dz = π · I_1(1) ≈ 1.7754662762...
    // (I_1 is the modified Bessel function of the first kind.)
    const double expected = M_PI * 0.5651591039924850;  // = I_1(1)

    qIntegral<LegendrePolynomial<32>> qint_leg;
    auto integrand_with_jac = [](const double& z) {
        return std::exp(z) * std::sqrt(1.0 - z*z);
    };
    const double via_legendre = qint_leg(integrand_with_jac, -1.0, 1.0);

    qIntegral<ChebyshevPolynomial2<32>> qint_cheb;
    auto integrand_no_jac = [](const double& z) { return std::exp(z); };
    const double via_chebyshev = qint_cheb(integrand_no_jac, -1.0, 1.0);

    // Chebyshev is essentially exact: it absorbs the √(1-z²) endpoint
    // singularity into the weights instead of multiplying it into the
    // integrand, so the remaining smooth factor exp(z) is integrated to
    // machine precision.
    CHECK(via_chebyshev == doctest::Approx(expected).epsilon(1e-12));

    // Legendre struggles with the √(1-z²) factor (endpoint derivative
    // singularity) — at order 32 the absolute error is ~5e-5. This is
    // exactly the motivation for the migration; the tolerance here is
    // loose to *document* the disparity, not to claim Legendre accuracy.
    CHECK(via_legendre  == doctest::Approx(expected).epsilon(1e-3));
    CHECK(via_chebyshev == doctest::Approx(via_legendre).epsilon(1e-3));
}

TEST_CASE("ChebyshevPolynomial2: 2D integration mirrors QPV/HVP usage pattern")
{
    // ∫_0^1 dx ∫_{-1}^{1} dz · x²·z²·√(1-z²) = (1/3) · (π/8) = π/24
    qIntegral2d<LegendrePolynomial<16>, ChebyshevPolynomial2<16>> qint2d;
    auto f = [](const double& x, const double& z) { return x*x * z*z; };
    const double res = qint2d(f, 0.0, 1.0, -1.0, 1.0);
    CHECK(res == doctest::Approx(M_PI / 24.0).epsilon(1e-12));
}
