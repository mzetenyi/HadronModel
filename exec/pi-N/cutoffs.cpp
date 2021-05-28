#include "cutoffs.hpp"
#include "units.hpp"
using namespace units_GeV;
#include "Config.hpp"

double cutoffRapp(halfint spin, double q2) {
  double Lambda(0.6*GeV);
  if (Config::exists("Lambda")) {
    Lambda = Config::get<double>("Lambda");
  }
  double L2(Lambda*Lambda);
  double exponent = 1.;
  if (spin==3*half) { exponent = 3./2.; }
  if (spin==5*half) { exponent = 5./2.; }
  return pow(L2/(L2+q2),exponent);
}

double cutoffZM(const string& resonance, double srt, double q2, double mr, double Gamma_Npi, int l) {
  static const double mn = Config::get<double>("Nucleon.mass");
  static const double mpi = Config::get<double>("pi_pm.mass");
  double mr2(mr*mr);
  double delta2;
  if (resonance == "D1232") {
    delta2 = POW<2>(0.3*GeV);
  } else if (resonance == "N1535") {
    delta2 = POW<2>(0.5*GeV);
  } else {
    delta2 = POW<2>(mr-mn-mpi) + POW<2>(Gamma_Npi)/4.;
  }
  double q02 = lambda(mr2,mn*mn,mpi*mpi)/(4.*mr2);
  return pow((q02+delta2)/(q2+delta2), (l+1.)/2.) * sqrt(mr/srt);
}

double cutoff_RNpi(const string& resonance, double srt) {
  static const double mn = Config::get<double>("Nucleon.mass");
  static const double mpi = Config::get<double>("pi_pm.mass");
  const double mr = Config::get<double>(resonance+".mass");
  const double Gamma_r = Config::get<double>(resonance+".width");
  const double BNpi = Config::get<double>(resonance+".BNpi");
  const int l = Config::get<int>(resonance+".l");
  const halfint spin = Config::get<halfint>(resonance+".spin");
  double Gamma_Npi = Gamma_r*BNpi;
  double q2 = lambda(srt*srt, mn*mn, mpi*mpi)/(4.*srt*srt);
  string cutType = Config::get<string>("cutoff");
  if (cutType == "Rapp") {
    return cutoffRapp(spin,q2);
  } else if (cutType == "ZM") {
    return cutoffZM(resonance,srt,q2,mr,Gamma_Npi,l);
  }
  cerr << "Unknown cutoff type in cutoff_RNpi: " << cutType << endl;
  exit(0);
}
