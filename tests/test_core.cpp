#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "Utils.hh"

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
