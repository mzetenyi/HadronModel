#ifndef CROSSECTION_HPP
#define CROSSECTION_HPP

#include "Config.hpp"
#include "Kinema.hpp"
#include "utils.hpp"

template <typename Amplitude>
class CrossSection {
 public:
  CrossSection(Amplitude amp) : amplitude(amp) {}
  double dsigma_dcosth(double costh) {
    Kinema2 kinIn = amplitude.getInputKinematics();
    Kinema2 kinOut = amplitude.getOutputKinematics();
    srt = kinIn.M;
    return 1. / (64. * pi_ * pi_ * srt * srt) * kinOut.pabs() / kinIn.pabs() *
           amplitude.MSQR();
  }
  double sigmaTot() {
    int ntheta = 20;
    if (isSet("CrossSection.ntheta"))
      ntheta = Config::get<int>("CrossSection.ntheta");
    double dtheta = 2./ntheta;
    double sum(0.);
    for (double theta(-1.+dtheta); theta<1.; theta+=dtheta) {
      sum += dsigma_dcosth(theta);
    }
    return sum/n * 2.;
  }

 private:
  Amplitude amplitude;
};

#endif  // CROSSECTION_HPP