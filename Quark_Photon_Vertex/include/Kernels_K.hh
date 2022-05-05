#pragma once
#include <cmath>
using namespace std;

double K11(double y,double l2, double u, double uprime, double V, double w, double wprime, double X)
{
	return -(1. + y * y) / 2 - y(1 - y)*X;
}



double K22(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
{
	return -(1. + y * y) / 2 (1.- 2. * l2 * V*V) + y(1. * y * y)* X ;
}



double K33(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
{
	return  y*(1. - 2. * l2 * V*V)  - (1.- y*y)* X;
}


double K44(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
{
	return y + (1. - y*y) * X ;
}


double K55(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
{
	return 3. * y ;
}


double K66(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
{
	return -y * (1. + 2. * l2 * V*V)*;
}


double K77(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
{
	return - y*y * (3. + 2. * l2 * V*V);
}

double K88(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
{
	return y * y - 2. * y * (1. - y*y) * X;
}

double K16(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
{
	return sqrt(2.) * (1. - y* y)* uprime * V  ;
}

double K61(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
{
	return - sqrt(2.) * (1. - y * y)* u * V;
}

double K17(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
{
	return -(1. - y*y)/sqrt(2.) * (1. + 2 * wprime * 2.*y * X)* uprime * V;
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
	return -(2. * y * uprime - (1 + y * y)* u) * V;;
}

double K67(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
{
	return 2. * y * ( uprime - y * u)*V ;
}

double K76(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
{
	return - 2. * y * (u - y * uprime)*V;

}

double K28(double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
{
	return K71(y, l2, u, uprime, V,  w,  wprime, X) 
		     + sqrt(2.)* (1. - y * y)  ;

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
	return  K55(y, l2, u, uprime, V, w, wprime, X)/y;

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




K(int i, int j, double y, double l2, double u, double uprime, double V, double w, double wprime, double X)
{
	if (i == 1 && j == 1)
	{
		return K11(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 2 && j == 2)
	{
		return K22(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 3 && j == 3)
	{
		return K33(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 4 && j == 4)
	{
		return K44(y, l2, u, uprime, V, w, wprime, X);
	}



	if (i == 5 && j == 5)
	{
		return K55(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 6 && j == 6)
	{
		return K66(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 7 && j == 7)
	{
		return K77(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 8 && j == 8)
	{
		return K88(y, l2, u, uprime, V, w, wprime, X);
	}


	if (i == 1 && j == 6)
	{
		return K16(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 6 && j == 1)
	{
		return K61(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 1 && j == 7)
	{
		return K17(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 7 && j == 1)
	{
		return K71(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 2 && j == 3)
	{
		return K23(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 3 && j == 2)
	{
		return K32(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 6 && j == 7)
	{
		return K67(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 7 && j == 6)
	{
		return K76(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 2 && j == 8)
	{
		return K28(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 8 && j == 2)
	{
		return K82(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 3 && j == 8)
	{
		return K38(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 8 && j == 3)
	{
		return K38(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 9 && j == 9)
	{
		return K99(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 10 && j == 10)
	{
		return K1010(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 10 && j == 11)
	{
		return K1011(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 11 && j == 10)
	{
		return K1110(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 11 && j == 11)
	{
		return K1111(y, l2, u, uprime, V, w, wprime, X);
	}

	if (i == 12 && j == 12)
	{
		return K1212(y, l2, u, uprime, V, w, wprime, X);

	}

	return 0;
}