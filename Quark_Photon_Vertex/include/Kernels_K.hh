#pragma once
#include <cmath>
using namespace std;

class K
{


	double K11(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return -(1. + y * y) / 2 - y * (1 - y)*X;
	}



	double K22(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return -(1. + y * y) / 2. *  (1. - 2. * l2 * V*V) + y(1. * y * y)* X;
	}



	double K33(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return  y * (1. - 2. * l2 * V*V) - (1. - y * y)* X;
	}


	double K44(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return y + (1. - y * y) * X;
	}


	double K55(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return 3. * y;
	}


	double K66(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return -y * (1. + 2. * l2 * V*V);
	}


	double K77(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return -y * y * (3. + 2. * l2 * V*V) + 2.*y*(1. - y * y) *X;
	}

	double K88(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return y * y - 2. * y * (1. - y * y) * X;
	}

	double K16(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return sqrt(2.) * (1. - y * y)* uprime * V;
	}

	double K61(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return -sqrt(2.) * (1. - y * y)* u * V;
	}

	double K17(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return -(1. - y * y) / sqrt(2.) * (1. + 2 * wprime * 2.*y * X)* uprime * V;
	}

	double K71(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return -(1. - y * y) / sqrt(2.) * (1. + 2 * w * 2.*y * X)* uprime * V;

	}

	double K23(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return (2. * y * u - (1 + y * y)* uprime) * V;
	}

	double K32(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return -(2. * y * uprime - (1. + y * y)* u) * V;
	}

	double K67(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return 2. * y * (uprime - y * u)*V;
	}

	double K76(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return -2. * y * (u - y * uprime)*V;

	}

	double K28(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return K71(y, l2, u, uprime, V, w, wprime, X)
			+ sqrt(2.)* (1. - y * y);

	}
	double K82(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return K17(y, l2, u, uprime, V, w, wprime, X)
			+ sqrt(2.)* (1. - y * y);

	}
	double K38(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return -K61(y, l2, u, uprime, V, w, wprime, X);

	}
	double K83(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return -K16(y, l2, u, uprime, V, w, wprime, X);

	}

	double K99(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return  K55(y, l2, u, uprime, V, w, wprime, X) / y;

	}


	double K1010(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return  K66(y, l2, u, uprime, V, w, wprime, X) / y;

	}

	double K1011(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return  K67(y, l2, u, uprime, V, w, wprime, X) / y;

	}

	double K1110(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return  K76(y, l2, u, uprime, V, w, wprime, X) / y;

	}

	double K1111(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return  K77(y, l2, u, uprime, V, w, wprime, X) / y;

	}

	double K1212(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
		return  K88(y, l2, u, uprime, V, w, wprime, X) / y;

	}

	using submatrixrow = vector<double>;
	using submatrix = vector<submatrixrow>;

	submatrix Kupper;
	submatrix Klower;




	public:
		K( double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
	{
			Kupper.resize(8, submatrixrow(8));

			Klower.resize(4, submatrixrow(4));


			Kupper[0][0] = K11(y, l2, u, uprime, V, w, wprime, X);

			Kupper[1][1] = K22(y, l2, u, uprime, V, w, wprime, X);

			Kupper[2][2] = K33(y, l2, u, uprime, V, w, wprime, X);

			Kupper[3][3] = K44(y, l2, u, uprime, V, w, wprime, X);

			Kupper[4][4] = K55(y, l2, u, uprime, V, w, wprime, X);

			Kupper[5][5] = K66(y, l2, u, uprime, V, w, wprime, X);

			Kupper[6][6] = K77(y, l2, u, uprime, V, w, wprime, X);

			Kupper[7][7] = K88(y, l2, u, uprime, V, w, wprime, X);

			Kupper[0][5] = K16(y, l2, u, uprime, V, w, wprime, X);

			Kupper[5][0] = K61(y, l2, u, uprime, V, w, wprime, X);

			Kupper[0][6] = K17(y, l2, u, uprime, V, w, wprime, X);
			Kupper[6][0] = K71(y, l2, u, uprime, V, w, wprime, X); 

			Kupper[5][6] = K67(y, l2, u, uprime, V, w, wprime, X);
			Kupper[6][5] = K76(y, l2, u, uprime, V, w, wprime, X);

			Kupper[1][2] = K23(y, l2, u, uprime, V, w, wprime, X);
			Kupper[2][1] = K32(y, l2, u, uprime, V, w, wprime, X);


			Kupper[1][7] = K28(y, l2, u, uprime, V, w, wprime, X);
			Kupper[7][1] = K82(y, l2, u, uprime, V, w, wprime, X);

			Kupper[2][7] = K38(y, l2, u, uprime, V, w, wprime, X);
			Kupper[7][2] = K83(y, l2, u, uprime, V, w, wprime, X);




			Klower[0][0]= K99(y, l2, u, uprime, V, w, wprime, X);
			Klower[1][1] = K1010(y, l2, u, uprime, V, w, wprime, X);
			Klower[1][2] = K1011(y, l2, u, uprime, V, w, wprime, X);
			Klower[2][1] = K1110(y, l2, u, uprime, V, w, wprime, X);
			Klower[2][2] = K1111(y, l2, u, uprime, V, w, wprime, X);
			Klower[3][3] = K1212(y, l2, u, uprime, V, w, wprime, X);

	}


		double get(int i, int j)
		{
			if (i < 8 && j < 8)
			{
				return Kupper[i][j];
			}

			else if (i >= 8 && j >= 8)
			{
				return Klower[i][j];
			}

			else {
				return 0;
			}
		}
}