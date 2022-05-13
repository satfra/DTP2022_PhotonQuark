#pragma once

#include <cmath>
#include <complex>

std::complex<double> g1(const std::complex<double> a10, const std::complex<double> a11, const double &z, const double &s)
{
    return a10 - (z / s) * a11;
}

std::complex<double> g2(const double &k, const std::complex<double> a11, const double &z, const double &s)
{
    return a11 / (k * k * z * s);
}

std::complex<double> g3(const double &k, const std::complex<double> a9,
          const double &z) //g3 it's actually immaginary and we shouldn't forget about that!!
{
    constexpr std::complex<double> II = {0.0, 1.0};
    return - II * a9 / (k * z);
}

std::complex<double> g4(const double &k, const std::complex<double> a12, const double &z, const double &s, const double &Q) //also g4 is immaginary
{
    constexpr std::complex<double> II = {0.0, 1.0};
    return II * a12 / (k * k * Q * s * z);
}

std::complex<double> f1(const double Q, const double s, const double z, const double a1, const double a6, const
double a7, const double a10,
   const double a11)
{
    return (a1 / (std::sqrt(2.) + z * s * (a6 + a11) - z * z * a7 - s * s * a10)) / (Q * Q * s * s);
}

std::complex<double> f2(const double Q, const double k, const double s, const double z, const double a2, const double
a3, const double a8,
   const double a12) // immaginari !!!!
{
    constexpr std::complex<double> II = {0.0, 1.0};
    return -II * (a2 / std::sqrt(2.) - a8 + s / z * (a3 / sqrt(2.) + a12)) / (k * k * Q * Q * Q * s * s);
}

std::complex<double> f3(const double Q, const double s, const double z, const double a2, const double a3) //immaginari
{
    constexpr std::complex<double> II = {0.0, 1.0};
    return II * (-a2 + z / s * a3) / (std::sqrt(2.) * Q);
}

std::complex<double> f4(const double Q, const double k, const double s, const double z, const double a4)
{
    return a4 / (std::sqrt(2.) * k * Q * s);
}

std::complex<double> f5(const double Q, const double k, const double s, const double z, const double a5, const double a9) //immaginari
{
    constexpr std::complex<double> II = {0.0, 1.0};
    return -II * (a5 - s / z * a9) / (k * s * Q * Q);
}

std::complex<double> f6(const double Q, const double k, const double s, const double z, const double a1, const double
a7, const double a11)
{
    return -(a1 / std::sqrt(2.) - a7 + s / z * a11) / (k * k * Q * Q * s * s);
}

std::complex<double> f7(const double Q, const double k, const double s, const double z, const double a1, const double
a6, const double a7)
{
    return -(a1 / std::sqrt(2.) - a7 + s / z * a6) / (k * k * Q * Q * s * s);
}

std::complex<double> f8(const double Q, const double k, const double s, const double a2, const double a8) // immaginario
{
    constexpr std::complex<double> II = {0.0, 1.0};
    return II * (a2 / std::sqrt(2.) - a8) / (k * k * Q * s * s);
}
