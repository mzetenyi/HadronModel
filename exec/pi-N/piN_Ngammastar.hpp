#ifndef PIN_NGAMMASTAR_HPP
#define PIN_NGAMMASTAR_HPP

#include <map>
#include <string>
#include "Vectors.hpp"
using namespace Vectors;
#include "Spinors.hpp"
using namespace Spinors;
#include "units.hpp"
using namespace units_GeV;
#include "Config.hpp"
#include "MultiArray.hpp"
#include "HelicityAmplitudes.hpp"
#include "Vrancx.hpp"
#include "observables.hpp"
#include "udouble.hpp"

class piN_Ngammastar {
public:
  piN_Ngammastar(double srt, double M);
  
  double u_Mandelstam(double costh) const;
  
  HelicityAmplitudes helicityAmplitudes(double costh) const;
  
  double cutoff_sch(const std::string& res) const;
  
  MultiArray<dcomplex> rho_prod(double costh) const;
  
  template <class Observable>
  MultiArray<double> partialExpectationValues(double costh, const Observable& OBS) const;
  
  template <class Observable>
  MultiArray<double> partialExpectationValues(const Observable& OBS) const;
  
  template <class Observable>
  udouble expectationValue(double costh, const Observable& OBS) const;

  template <class Observable>
  udouble expectationValue(const Observable& OBS) const;

  vector<string> getResonances() const;
  vector<string> getChannels() const;
  vector<string> getEMchannels() const;
private:
  const double srt; ///< srt(s)
  const double M; ///< dilepton mass
  vector<string> resonances;
  vector<string> channels; // Born,resonances
  vector<string> EMchannels; // gamma, rho
  vector<string> DECchannels; // pi, gamma, rho
  NamedIndex NIch;
  NamedIndex NIem;
  NamedIndex NIdec;
  //  vector<NamedIndex> ChI; ///< channel index by name
  //vector<string> couplings;
  //NamedIndex couplIDX;
  MultiArray<double> relative_errors;
  // switches of different channels:
  bool gamma;
  bool rho;
  bool rho2;
  bool VMD2; // if true: use VMD2 version (in contrast to Kroll-Lee-Zumino)
  bool vRho;
  bool tRho;
  bool vGamma;
  bool tGamma;
  bool vOmega;
  bool tOmega;
  bool sch;
  bool uch;
  bool Born;
  bool Born_s;
  bool Born_u;
  bool Born_t;
  bool Born_c;
  bool Born_corr;
  double s;
  double M2;
  double mpi_pm;
  double mpi_pm2;
  double mpi_0;
  double mpi;
  double mpi2;
  double mrho;
  double mn;
  double mn2;
  double EN_in;
  double Epi_in;
  double pin_abs;
  FourVector p1;
  FourVector p2;
  FourVector q;
  FourVector P;
  double EN_out;
  double k0;
  double pout_abs;
  dcomplex F_rho;
  dcomplex F_ome;
  vector<string> determineResonances() const;
  vector<string> determineChannels() const;
  vector<string> determineEMchannels() const;
  // vector<NamedIndex> getChannelIndices() const;
  // vector<string> getCouplings() const;
  MultiArray<double> getRelative_errors() const;
};

template <class Observable>
MultiArray<double> piN_Ngammastar::partialExpectationValues(double costh, const Observable& OBS) const {
  int npol = 2;
  double fac = 1./(64.*POW<4>(2.*pi_)*s) * pout_abs/pin_abs * M * 1./npol;
  MultiArray<double> PartExpVal({NIch,NIem});
  HelicityAmplitudes HA = helicityAmplitudes(costh);
  for (auto R : channels) {
    for (auto al : EMchannels) {
      dcomplex O_Ral(0);
      for (auto Rp : channels) {
        for (auto alp : EMchannels) {
          for (halfint la : {-_1,_0,_1}) {
            for (halfint lap : {-_1,_0,_1}) {
              dcomplex rho_la_lap(0);
              for (halfint la1 : {-half,half}) {
                for (halfint la2 : {-half,half}) {
                  rho_la_lap += HA({R,al},{la1,la2,la}) * conj(HA({Rp,alp},{la1,la2,lap}))
                    + HA({Rp,alp},{la1,la2,la}) * conj(HA({R,al},{la1,la2,lap}));
                }
              }
              O_Ral += 1./2. * fac * rho_la_lap * OBS(lap,la);
            }
          }
        }
      }
      PartExpVal(NIch(R),NIem(al)) = real(O_Ral);
    }
  }
  return PartExpVal;
}


template <class Observable>
MultiArray<double> piN_Ngammastar::partialExpectationValues(const Observable& OBS) const {
  int npol = 2;
  double fac = 1./(64.*POW<4>(2.*pi_)*s) * pout_abs/pin_abs * M * 1./npol;
  //  cerr << "fac = " << fac << endl;
  MultiArray<double> PartExpVal({NIch,NIem});

  uint nth(20);
  if (Config::exists("nth")) nth = Config::get<int>("nth");
  const double startcosth(-1.);
  const double maxcosth(1.);
  const double dcosth = (maxcosth-startcosth)/nth;
  
  for (auto R : channels) {
    for (auto al : EMchannels) {
      PartExpVal(NIch(R),NIem(al)) = 0;
    }
  }
  
  for (double costh(startcosth+dcosth/2.); costh<maxcosth; costh+=dcosth) {
    HelicityAmplitudes HA = helicityAmplitudes(costh);
    for (auto R : channels) {
      for (auto al : EMchannels) {
        dcomplex O_Ral(0);
        //        dcomplex O_Ral_test(0);
        for (auto Rp : channels) {
          for (auto alp : EMchannels) {
            for (halfint la : {-_1,_0,_1}) {
              for (halfint lap : {-_1,_0,_1}) {
                dcomplex rho_la_lap(0);
                for (halfint la1 : {-half,half}) {
                  for (halfint la2 : {-half,half}) {
                    rho_la_lap += HA({R,al},{la1,la2,la}) * conj(HA({Rp,alp},{la1,la2,lap}))
                      + HA({Rp,alp},{la1,la2,la}) * conj(HA({R,al},{la1,la2,lap}));
                  }
                }
                O_Ral += 1./2. * fac * rho_la_lap * OBS(lap,la);
                //      O_Ral_test += 1./2. * rho_la_lap * OBS(lap,la);
                //cerr << "O_Ral/O_Ral_test = " << O_Ral/O_Ral_test << endl;
              }
            }
          }
        }
        PartExpVal(NIch(R),NIem(al)) += dcosth * real(O_Ral);
      }
    }
  }
  // for (auto R : channels) {
  //   for (auto al : EMchannels) {
  //     cerr << R << al << " -> " << PartExpVal(NIch(R),NIem(al)) << endl;
  //   }
  // }
  return PartExpVal;
}


template <class Observable>
udouble piN_Ngammastar::expectationValue(double costh, const Observable& OBS) const {
  MultiArray<double> PartExpVal = partialExpectationValues(costh,OBS);
  // for (auto R : channels) {
  //   for (auto al : EMchannels) {
  //     cerr << R << al << " -> " << PartExpVal(NIch(R),NIem(al)) << endl;
  //   }
  // }
  double expVal(0);
  double r2(0);
  for (auto R : channels) {
    double rpi = relative_errors(NIch(R),NIdec("pi"));
    double sum_O2(0);
    double sum_r2O2(0);
    for (auto al : EMchannels) {
      double r = relative_errors(NIch(R),NIdec(al));
      double O_Ral = PartExpVal(NIch(R),NIem(al));
      expVal += O_Ral;
      sum_O2 += O_Ral*O_Ral;
      sum_r2O2 += r*r * O_Ral*O_Ral;
    }
    r2 += rpi*rpi * sum_O2 + sum_r2O2;
  }
  r2 *= 4./(expVal*expVal);
  return udouble(expVal,sqrt(r2));
}
  
  
template <class Observable>
udouble piN_Ngammastar::expectationValue(const Observable& OBS) const {
  MultiArray<double> PartExpVal = partialExpectationValues(OBS);
  // for (auto R : channels) {
  //   for (auto al : EMchannels) {
  //     cerr << R << al << " -> " << PartExpVal(NIch(R),NIem(al)) << endl;
  //   }
  // }
  double expVal(0);
  double r2(0);
  for (auto R : channels) {
    double rpi = relative_errors(NIch(R),NIdec("pi"));
    double sum_O2(0);
    double sum_r2O2(0);
    for (auto al : EMchannels) {
      double r = relative_errors(NIch(R),NIdec(al));
      double O_Ral = PartExpVal(NIch(R),NIem(al));
      expVal += O_Ral;
      sum_O2 += O_Ral*O_Ral;
      sum_r2O2 += r*r * O_Ral*O_Ral;
    }
    r2 += rpi*rpi * sum_O2 + sum_r2O2;
  }
  r2 *= 4./(expVal*expVal);
  //  cerr << "return from expectationValue: " << expVal << " (" << sqrt(r2) << ")" << endl; 
  return udouble(expVal,sqrt(r2));
}
  
  
#endif // PIN_NGAMMASTAR_HPP
