#pragma once //include stuff only once 

#include <cmath>
#include <Utils.hh>
#include <vector>
#include <complex>

static constexpr std::complex<double> II = { 0.0, 1.0 };

namespace basistransform
{
	std::complex<double> f1(const double Q, const double s, const double z,const double k, std::vector< std::complex<double> > a)
	{
		return 1.0 / powr<2>(Q*s) * (a[0] /std::sqrt(2.) + z * s * (a[5] + a[10]) - powr<2>(z) * a[6] - powr<2>(s) * a[9]);
	}
	
	std::complex<double> f2(const double Q, const double s, const double z,const double k, std::vector< std::complex<double> > a)
	{
		return II / (powr<3>(Q)*powr<2>(k*s)) * (a[1] / std::sqrt(2.) - a[7] + s/z * (a[2]/std::sqrt(2.) + a[11]));
	}


	std::complex<double> f3(const double Q, const double s, const double z,const double k, std::vector< std::complex<double> > a)
	{
		return - II / (std::sqrt(2.) * Q) * (-a[1] + z/s * a[2] );
	}

	std::complex<double> f4(const double Q, const double s, const double z, const double k, std::vector<std::complex<double>> a)
	{
		return 1.0 / (std::sqrt(2.) * k * s * Q) * a[3];
	}

	std::complex<double> f5(const double Q, const double s, const double z,const double k, std::vector< std::complex<double> > a)
	{
		return -II / (powr<2>(Q)*k*s) * (a[4]  - s/z * a[8]);
	}

	std::complex<double> f6(const double Q, const double s, const double z,const double k, std::vector<std::complex<double> > a)
	{
		return -1.0/powr<2>(Q*s*k) * ( a[0]/std::sqrt(2) - a[6] + s/z * a[10]);
	}

	std::complex<double> f7(const double Q, const double s, const double z, const double k,std::vector< std::complex<double> > a)
	{
		return -1.0 / powr<2>(Q*s*k) * (a[0] / sqrt(2.) - a[6] + s/z * a[5]);
	}

	std::complex<double> f8(const double Q, const double s, const double z, const double k,std::vector< std::complex<double> > a)
	{
		return - II  / ( Q *powr<2>(k*s) ) * (a[1] / sqrt(2.) - a[7]);
	}


	std::complex<double> g1(const double Q, const double s, const double z,const double k, std::vector< std::complex<double> > a)
	{
		return a[9]- z/s * a[10] ;
	}

	std::complex<double> g2(const double Q, const double s, const double z,const double k, std::vector< std::complex<double> > a)
	{
		return 1. / (powr<2>(k)*z*s) * a[10];
	}

	std::complex<double> g3(const double Q, const double s, const double z, const double k,std::vector< std::complex<double> > a)
	{
		return II/(k*z) * a[8];
	}

	std::complex<double> g4(const double Q, const double s, const double z, const double k,std::vector< std::complex<double> > a)
	{
		return II/(powr<2>(k)*z *s) * a[11];
	}
}
