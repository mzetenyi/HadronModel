#ifndef AMPLITUDE_PIN_ELASTIC_HPP
#define AMPLITUDE_PIN_ELASTIC_HPP

#include "Kinema.hpp"

class Amplitude_piN_elastic {
 public:
  Amplitude_piN_elastic(double srt, int QN_in, int Qpi_in, int QN_out,
                        int Qpi_out);
  Kinema2 getInputKinematics();
  Kinema2 getOutputKinematics();
  double MSQR(double costh);

 private:
  Amplitude_piN_elastic();
  const double srt;
  const int QN_in;
  const int Qpi_in;
  const int QN_out;
  const int Qpi_out;
  const double mN;
  const double mpi_in;
  const double mpi_out;
  const Kinema2 kinIn;
  const Kinema2 kinOut;
};

#endif  // AMPLITUDE_PIN_ELASTIC_HPP