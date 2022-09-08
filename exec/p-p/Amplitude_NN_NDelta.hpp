#ifndef NN_NDELTA_HPP
#define NN_NDELTA_HPP

#include "Kinema.hpp"

class Amplitude_NN_NDelta {
 public:
  Amplitude_NN_NDelta(double srt, int QN, int QD);
  Kinema2 getInputKinematics();
  Kinema2 getOutputKinematics(double m);
  double width(double m);
  double spectralFunction(double m);
  double MSQR(double m, double costh);

 private:
  Amplitude_NN_NDelta();
  double srt;
  int QN;
  int QD;
};

#endif  // NN_NDELTA_HPP
