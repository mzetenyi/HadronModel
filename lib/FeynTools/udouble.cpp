#include "udouble.hpp"
#include <iostream>
#include <iomanip>

using namespace std;

udouble::udouble(double value, double rel_uncert) :
  value(value),
  rel_uncert(rel_uncert) {
  }

double udouble::get_value() const {
  return value;
}

double udouble::get_rel_uncert() const {
  return rel_uncert;
}

double udouble::get_uncert() const {
  return rel_uncert * value;
}

double udouble::minimum() const {
  return value*(1. - rel_uncert);
}

double udouble::maximum() const {
  return value*(1. + rel_uncert);
}

udouble operator*(double x, const udouble& u) {
  return udouble(x*u.value,u.rel_uncert);
}

ostream& operator<<(ostream& out, const udouble& ud) {
  out << setw(14) << ud.value
      << setw(14) << ud.rel_uncert 
      << setw(14) << ud.minimum()
      << setw(14) << ud.maximum();
  return out;
}

