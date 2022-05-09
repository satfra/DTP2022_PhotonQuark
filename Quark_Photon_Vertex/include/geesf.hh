#pragma once
#include <cmath>

// your code
using namespace std

double g1 (const double& a10, const double& a11, const double& z, const double& s)
{
 return a10-(z/s)*a11;
}

double g2 (const double& k, const double& a11, const double& z, const double& s)
{
 return a11/(k*k*z*s);
}

double g3 (const double& k, const double& a9, const double& z) //g3 it's actually immaginary and we shouldn't forget about that!!
{
 return a9/(k*z);
}

double g4 (const double& k, const double& a12, const double& z, const double& s, const double& Q) //also g4 is immaginary
{
 return a12/(k*k*Q*s*z);
}


double f1 (const double Q, const double s, const double z, const double a1, const double a6, const double a7 , const double a10, const double a11)
{
	return (a1 / (sqrt(2) + z * s * (a6 + a11) - z * z * a7 - s * s * a10) / (Q * Q * s * s);
} 

double f2 (const double Q, const double k, const double s, const double z, const double a2, const double a3, const double a8 , const double a12) // immaginari !!!!
{
	return -(a2 / (sqrt(2) - a8 + s / z * (a3/ (sqrt(2)) + a12 ) / (k * k * Q* Q * Q * s * s);
}

double f3 (const double Q, const double s, const double z, const double a2, const double a3) //immaginari
{
	return (-a2 + z / s * a3) / (sqrt(2) * Q);
}

double f4 (const double Q, const double k, const double s, const double z, const double a4)
{
	return a4 / (sqrt(2) * Q * k * s);
}

double f5 (const double Q, const double k, const double s, const double z, const double a5, const double a9) //immaginari
{
	return (a5 - s / z * a9) / (k * s * Q * Q);
}

double f6 (const double Q, const double k, const double s, const double z, const double a1, const double a7, const double a11)
{
	return -(a1 / (sqrt(2)) - a7 + s / z * a11) / (k * k * Q * Q * s * s);
}

double f7 (const double Q, const double k, const double s, const double z, const double a1, const double a6, const double a7)
{
	return -(a1 / (sqrt(2)) - a7 + s / z * a6) / (k * k * Q * Q * s * s);
}

double f8 (const double Q, const double k, const double s, const double a2, const double a8) // immaginario
{
	return (a2 / (sqrt(2)) - a8) / (k * k * Q * s * s);
}
