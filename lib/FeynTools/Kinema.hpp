#ifndef KINEMA_HPP
#define KINEMA_HPP

#include "Vectors.hpp"
using namespace Vectors;

double momentum(double M, double m1, double m2);

class Kinema2 {
 public:
  Kinema2(double M, double m1, double m2);
  double getm(int i) const;
  double getpabs() const;
  double E1() const;
  double E2() const;
  double getE(int i) const;
  double threshold() const;
  FourVector P() const;
  FourVector p1(double costh = 0, double phi = 0) const;
  FourVector p2(double costh = 0, double phi = 0) const;
  FourVector getp(int i, double costh = 0, double phi = 0) const;
  FourVector transformFrom(int i, double costh, double phi, const FourVector& q) const; 
  const double M;
  const double m1;
  const double m2;
};

#endif  // KINEMA_HPP
