#include "Kinema.hpp"

#include <cassert>

#include "utils.hpp"

double momentum(double M, double m1, double m2) {
  if (M < m1 + m2) return 0;
  return sqrt(lambda(M * M, m1 * m1, m2 * m2)) / (2. * M);
}

Kinema2::Kinema2(double M, double m1, double m2) : M(M), m1(m1), m2(m2) {}

double Kinema2::getm(int i) const {
  assert(i == 1 or i == 2);
  return (i == 1) ? m1 : m2;
}

double Kinema2::getpabs() const {
  if (M < m1 + m2) return 0;
  return sqrt(lambda(M * M, m1 * m1, m2 * m2)) / (2. * M);
}

double Kinema2::E1() const { return (M * M + m1 * m1 - m2 * m2) / (2. * M); }

double Kinema2::E2() const { return (M * M - m1 * m1 + m2 * m2) / (2. * M); }

double Kinema2::getE(int i) const {
  assert(i == 1 or i == 2);
  return (i == 1) ? E1() : E2();
}

double Kinema2::threshold() const { return m1 + m2; }

FourVector Kinema2::P() const { return FourVector(M, 0, 0, 0); }

FourVector Kinema2::p1(double costh, double phi) const {
  double pp = getpabs();
  double sinth = sqrt(1. - costh * costh);
  return FourVector(E1(), pp * sinth * cos(phi), pp * sinth * sin(phi),
                    pp * costh);
}

FourVector Kinema2::p2(double costh, double phi) const {
  double pp = getpabs();
  double sinth = sqrt(1. - costh * costh);
  return FourVector(E2(), -pp * sinth * cos(phi), -pp * sinth * sin(phi),
                    -pp * costh);
}

FourVector Kinema2::getp(int i, double costh, double phi) const {
  assert(i == 1 or i == 2);
  return (i == 1) ? p1(costh, phi) : p2(costh, phi);
}

FourVector Kinema2::transformFrom(int i, double costh, double phi,
                                  const FourVector& q) const {
  assert(i == 1 or i == 2);
  double pabs = getpabs();
  double E = getE(i);
  double m = getm(i);
  double ch = E / m;
  double sh = pabs / m;
  double costh_ = (i == 1) ? costh : -costh;
  double sinth_ = sqrt(1. - costh_ * costh_);
  double phi_ = (i == 1) ? phi : phi + pi_;
  FourVector q1(ch * q(0) + sh * q(3), q(1), q(2), sh * q(0) + ch * q(3));
  FourVector q2(q1(0), q1(1) * costh_ + q1(3) * sinth_, q1(2),
                - q1(1) * sinth_ + q1(3) * costh_);
  FourVector q3(q2(0), q2(1) * cos(phi_) - q2(2) * sin(phi_),
                q2(1) * sin(phi_) + q2(2) * cos(phi_), q2(3));
  return q3;
}