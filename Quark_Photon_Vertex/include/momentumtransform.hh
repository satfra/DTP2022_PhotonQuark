#pragma once //include stuff only once 
#include <math.h>
#include <Utils.hh>

double l2 (double k, double k1, double z, double z1, double y)
{
 return k*k+k1*k1-2*k*k1*(z*z1+y*sqrt(1-z*z)*sqrt(1-z1*z1));
} 

double u (double k, double z)
{
 return k * sqrt (1 - powr<2>(z));
}

double V (double k, double z, double k1, double z1, double l)
{
 return (k * z - k1 * z1) / (powr<2>(l))
}

double w (double u, double l)
{
 return powr<2>(u/l);
}

double X (double u, double u1, double l)
{
 return u*u1/powr<2>(l);
}


int main(int argc, char* argv[])
{
  int i;
  double u=1, u'=2, z=2, k, k'; 
 
 
   printf("for u=%f u'= %f l=%f we get X=%f", u, u', z, X(u, u', z);
  
  
  
  return 0;
}
