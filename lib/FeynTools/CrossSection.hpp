#ifndef CROSSECTION_HPP
#define CROSSECTION_HPP

#include "Config.hpp"
#include "Kinema.hpp"
#include "utils.hpp"

template <typename Amplitude>
class CrossSectionSS {
 public:
  CrossSectionSS(Amplitude amp) : amplitude(amp) {}
  double dsigma_dcosth(double costh) {
    Kinema2 kinIn = amplitude.getInputKinematics();
    Kinema2 kinOut = amplitude.getOutputKinematics();
    double srt = kinIn.M;
    return 1. / (64. * pi_ * pi_ * srt * srt) * kinOut.pabs() / kinIn.pabs() *
           amplitude.MSQR();
  }
  double sigmaTot() {
    int ntheta = 20;
    if (isSet("CrossSection.ntheta"))
      ntheta = Config::get<int>("CrossSection.ntheta");
    double dtheta = 2. / ntheta;
    double sum(0.);
    for (double theta(-1. + dtheta); theta < 1.; theta += dtheta) {
      sum += dsigma_dcosth(theta);
    }
    return sum / ntheta * 2.;
  }

 private:
  Amplitude amplitude;
};

template <typename Amplitude>
class CrossSectionSR {
 public:
  CrossSectionSR(Amplitude amp) : amplitude(amp) {}
  double dsigma_dm_dcosth(double m, double costh) {
    Kinema2 kinIn = amplitude.getInputKinematics();
    Kinema2 kinOut = amplitude.getOutputKinematics(m);
    double srt = kinIn.M;
    return 1. / (64. * pi_ * pi_ * srt * srt) * kinOut.pabs() / kinIn.pabs() *
           amplitude.MSQR(m, costh) * 2. * m * amplitude.spectralFunction(m);
  }
  double dsigma_dm(double m) {
    int ncosth = 20;
    if (isSet("CrossSection.ntheta"))
      ncosth = Config::get<int>("CrossSection.ncosth");
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
    int nm = 20;
    if (isSet("CrossSection.nm")) nm = Config::get<int>("CrossSection.nm");
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