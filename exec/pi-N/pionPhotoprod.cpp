#include "pionPhotoprod.hpp"
#include "piN_Ngammastar.hpp"

udouble pionPhotoprod_sigtot(double srt) {
  piN_Ngammastar DS(srt,0.000001*MeV);
  Observable_Xsec_PionPhotoprod OBS;
  udouble ret = DS.inverseExpectationValue(OBS);
  return ret;
}