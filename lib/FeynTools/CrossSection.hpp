#ifndef CROSSECTION_HPP
#define CROSSECTION_HPP

#include "Config.hpp"
#include "Kinema.hpp"
#include "utils.hpp"

/**
 * @brief 2 -> 2 cros section with two stable particles in the final state.
 * 
 * @tparam Amplitude 
 */
template <typename Amplitude>
class CrossSectionSS {
 public:
  CrossSectionSS(Amplitude amp) : amplitude(amp) {}
  double dsigma_dcosth(double costh) {
    Kinema2 kinIn = amplitude.getInputKinematics();
    Kinema2 kinOut = amplitude.getOutputKinematics();
    double srt = kinIn.M;
    return 1. / (32. * pi_ * srt * srt) * kinOut.getpabs() / kinIn.getpabs() *
           amplitude.MSQR(costh);
  }
  double sigmaTot() {
    int ncosth = getParam<int>("ncosth",20);
    double dcosth = 2. / ncosth;
    double sum(0.);
    for (double costh(-1. + dcosth); costh < 1.; costh += dcosth) {
      sum += dsigma_dcosth(costh);
    }
    return sum / ncosth * 2.;
  }

 private:
  Amplitude amplitude;
};

/**
 * @brief 2 -> 2 cross section with a stable particle and a resonance in the final state.
 * 
 * @tparam Amplitude 
 */
template <typename Amplitude>
class CrossSectionSR {
 public:
  CrossSectionSR(Amplitude amp) : amplitude(amp) {}
  double dsigma_dm_dcosth(double m, double costh) {
    Kinema2 kinIn = amplitude.getInputKinematics();
    Kinema2 kinOut = amplitude.getOutputKinematics(m);
    double srt = kinIn.M;
    return 1. / (32. * pi_ * srt * srt) * kinOut.getpabs() / kinIn.getpabs() *
           amplitude.MSQR(m, costh) * 2. * m * amplitude.spectralFunction(m);
  }
  double dsigma_dm(double m) {
    int ncosth = getParam<int>("ncosth",20);
    double dcosth = 2. / ncosth;
    double sum(0.);
    for (double costh(-1. + dcosth); costh < 1.; costh += dcosth) {
      sum += dsigma_dm_dcosth(m, costh);
    }
    return sum / ncosth * 2.;
  }
  double sigmaTot() {
    Kinema2 kinOut = amplitude.getOutputKinematics(0.);
    double srt = kinOut.M;
    double m1 = kinOut.m1;
    double mmax = srt - m1;
    int nm = getParam<int>("nm",20);
    double dm = 2. / nm;
    double sum(0.);
    for (double m(dm); m < mmax; m += dm) {
      sum += dsigma_dm(m);
    }
    return sum / nm * mmax;
  }

 private:
  Amplitude amplitude;
};

#endif  // CROSSECTION_HPP