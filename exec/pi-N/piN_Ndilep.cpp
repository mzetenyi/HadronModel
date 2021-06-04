#include "utils.hpp"
#include "udouble.hpp"
#include "piN_Ndilep.hpp"
#include "piN_Ngammastar.hpp"

udouble piN_Ndilep_dsigma_dM(double srt, double M) {
  piN_Ngammastar DS(srt,M);
  Observable_Xsec_Dilepton_dM OBS(M);
  udouble ret = M/pi_*DS.expectationValue(OBS);
  return ret;
}

udouble piN_Ndilep_dsigma_dM_dcosth(double srt, double M, double costh) {
  piN_Ngammastar DS(srt,M);
  Observable_Xsec_Dilepton_dM OBS(M);
  udouble ret = M/pi_*DS.expectationValue(costh,OBS);
  return ret;
}

double piN_Ndilep_sigma_highM(double srt, double Mmin) {
  static const double mn = Config::get<double>("Nucleon.mass");
  double Mmax(srt-mn);

  int NM(100);
  if (Config::exists("NM")) {
    NM = Config::get<int>("NM");
  }
  double dM = (Mmax-Mmin)/NM;
  double sum(0.);
  for (double M(Mmin+dM/2.); M<Mmax; M+=dM) {
    sum += piN_Ndilep_dsigma_dM(srt,M).get_value();
  }
  return sum*dM;
}
