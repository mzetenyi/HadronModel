#ifndef AMPLITUDE_R_PIN_HPP
#define AMPLITUDE_R_PIN_HPP

#include <string>

#include "Kinema.hpp"

class Amplitude_R_piN {
 public:
  Amplitude_R_piN(std::string resonance, double srt, int QR, int QN, int Qpi);
  double MSQR();
  Kinema2 getKinematics();

 private:
  Amplitude_R_piN();

 protected:
  std::string resonance;

 public:
  const double M;

 protected:
  double srt;
  int QR;
  int QN;
  int Qpi;
  double mN;
  double mpi;
  Kinema2 kin;
};

#endif  // AMPLITUDE_R_PIN_HPP