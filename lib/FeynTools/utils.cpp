#include "utils.hpp"

bool isnan(const dcomplex& z) { return (std::isnan(z.real()) or std::isnan(z.imag())); } 

double isreal(const dcomplex& z, double eps) {
  //  if (z.real() == 0) { return z.imag() == 0; }
  return fabs(z.imag()) < eps;
}

template <>
double POW<0>(double x) {
  return 1;
}

template <>
dcomplex POW<0>(dcomplex x) {
  return 1;
}

double lambda(double x, double y, double z) {
  return x*x + y*y + z*z - 2.*(x*y + y*z + z*x);
}

double myAcos(double x) {
  return (x>1.) ? 0 : ((x<-1.) ? pi_ : acos(x));
 }

double myAcos(double x, double y) {
  return (y==0.) ? 0 : myAcos(x/y);
}
