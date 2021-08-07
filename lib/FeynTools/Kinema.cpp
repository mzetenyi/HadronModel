#include "Kinema.hpp"

#include "utils.hpp"

Kinema2::Kinema2(double M, double m1, double m2) : M(M), m1(m1), m2(m2) {}

double Kinema2::pabs() const {
  return sqrt(lambda(M * M, m1 * m1, m2 * m2)) / (2. * M);
}

double Kinema2::E1() const { return (M * M + m1 * m1 - m2 * m2) / (2. * M); }

double Kinema2::E2() const { return (M * M - m1 * m1 + m2 * m2) / (2. * M); }

FourVector Kinema2::p1(double costh, double phi) const {
  double pp = pabs();
  double sinth = sqrt(1. - costh*costh);
  return FourVector(E1(), pp * costh * cos(phi),
                    pp * costh * sin(phi), pp * sinth);
}

FourVector Kinema2::p2(double costh, double phi) const {
  double pp = pabs();
  double sinth = sqrt(1. - costh*costh);
  return FourVector(E1(), -pp * costh * cos(phi),
                    -pp * costh * sin(phi), -pp * sinth);
}