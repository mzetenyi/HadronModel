#include "Transformations.hpp"

#include <iostream>

#include "utils.hpp"

FourTensor LorentzBoostX(double khi) {
  return FourTensor::createUpDown(cosh(khi), sinh(khi), 0, 0, sinh(khi),
                                  cosh(khi), 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
}

FourTensor LorentzBoostY(double khi) {
  return FourTensor::createUpDown(cosh(khi), 0, sinh(khi), 0, 0, 1, 0, 0,
                                  sinh(khi), 0, cosh(khi), 0, 0, 0, 0, 1);
}

FourTensor LorentzBoostZ(double khi) {
  return FourTensor::createUpDown(cosh(khi), 0, 0, sinh(khi), 0, 1, 0, 0, 0, 0,
                                  1, 0, sinh(khi), 0, 0, cosh(khi));
}

FourTensor RotationX(double theta) {
  return FourTensor::createUpDown(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, cos(theta),
                                  -sin(theta), 0, 0, sin(theta), cos(theta));
}

FourTensor RotationY(double theta) {
  return FourTensor::createUpDown(1, 0, 0, 0, 0, cos(theta), 0, sin(theta), 0,
                                  0, 1, 0, 0, -sin(theta), 0, cos(theta));
}

FourTensor RotationZ(double theta) {
  return FourTensor::createUpDown(1, 0, 0, 0, 0, cos(theta), -sin(theta), 0, 0,
                                  sin(theta), cos(theta), 0, 0, 0, 0, 1);
}

FourTensor TransformationTo(FourVector p) {
  if (not p.timelike()) {
    std::cerr << "TransformationTo(FourVector p) - p is not timelike: p = " << p
              << std::endl;
    exit(0);
  }
  if (p.atRest()) {
    return FourTensor::unitTensor;
  }
  double m = sqrt(p * p);
  double khi = cosh(p(0) / m);
  ThreeVector P = p.spacial();
  double Pabs = sqrt(P * P);
  double theta = myAcos(P(3) / Pabs);
  double P12abs = Pabs * sin(theta);
  double phi = myAcos(P(1) / P12abs);
  return LorentzBoostZ(-khi) * RotationY(-theta) * RotationZ(-phi);
}

FourTensor TransformationFrom(FourVector p) {
  if (not p.timelike()) {
    std::cerr << "TransformationFrom(FourVector p) - p is not timelike: p = "
              << p << std::endl;
    exit(0);
  }
  if (p.atRest()) {
    return FourTensor::unitTensor;
  }
  double m = sqrt(p * p);
  double khi = acosh(p(0) / m);
  ThreeVector P = p.spacial();
  double Pabs = sqrt(P * P);
  double theta = myAcos(P(3), Pabs);
  double P12abs = Pabs * sin(theta);
  double phi = myAcos(P(1), P12abs);
  return RotationZ(phi) * RotationY(theta) * LorentzBoostZ(khi);
}
