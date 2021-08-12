#include "pionPhotoprodAnalytic.hpp"

#include "Config.hpp"
#include "Gamma.hpp"
#include "Kinema.hpp"
#include "Vectors.hpp"

#include <iostream>

using namespace std;

pionPhotoprodAnalytic::pionPhotoprodAnalytic(double srt)
    : srt(srt),
      mN(Config::get<double>("Nucleon.mass")),
      mpi(Config::get<double>("pi_pm.mass")),
      KINin(srt, mN, 0),
      KINout(srt, mN, mpi) {}

pionPhotoprodAnalytic::~pionPhotoprodAnalytic() {}

double pionPhotoprodAnalytic::MSQR(double costh) {
  double mR = Config::get<double>("N1440.mass");
  FourVector pi = KINin.p1();
  FourVector k = KINin.p2();
  FourVector pf = KINout.p1();
  double pi_k = pi * k;
  double pf_k = pf * k;
  double pi_pf = pi * pf;
  double mN2 = mN * mN;
  double mN3 = mN2 * mN;
  double mN4 = mN3 * mN;
  double mR2 = mR * mR;
  double s = srt * srt;

  int npol = 2 * 2;  // initial pol. (photon, nucleon)
  double iso = sqrt(2.);
  double gRNg = Config::get<double>("N1440.gngamma");
  double gRNpi = Config::get<double>("N1440.g0");
  double Gr = Config::get<double>("N1440.width");
  double mrho = Config::get<double>("rho.mass");
  const int l = Config::get<int>("N1440.l");
  const dcomplex BWs = BW(s, mR, Gamma_R("N1440", s, mR, Gr, l));
  double BWs2 = POW<2>(abs(BWs));
  double fac = iso * gRNpi * gRNg / (2. * mpi * mrho);

  double ret =
      128. * POW<2>(fac) * 1. / npol * BWs2 *
      (+pi_pf * POW<2>(pi_k) * mN2 + 2. * pi_pf * POW<2>(pi_k) * mR * mN +
       pi_pf * POW<2>(pi_k) * mR2 + 2. * POW<2>(pi_k) * pf_k * mN2 +
       2. * POW<2>(pi_k) * pf_k * mR * mN - POW<2>(pi_k) * mN4 -
       2. * POW<2>(pi_k) * mR * mN3 - POW<2>(pi_k) * mR2 * mN2 +
       2. * POW<3>(pi_k) * pf_k - 2. * POW<3>(pi_k) * mN2 -
       2. * POW<3>(pi_k) * mR * mN);
  cerr << "costh = " << costh << endl;
  PR(mR); PR(gRNpi); PR(gRNg);
  cerr << "fac = " << fac << endl;
  //cerr << "BWs = " << BWs << endl;
  cerr << "MSQR = " << ret << endl;
  return ret;
}

double pionPhotoprodAnalytic::dsig_dcosth(double costh) {
  double fac = 1. / (32. * pi_ * srt * srt) * KINout.pabs() / KINin.pabs();
  cerr << "s = " << srt*srt << endl;
  cerr << "pout = " << KINout.pabs() << endl;
  cerr << "pin  = " << KINin.pabs() << endl;
  cerr << "kine fac = " << fac << endl;
  return fac * MSQR(costh);
};

double pionPhotoprodAnalytic::sigtot() {
  double dcosth = 0.01;
  if (Config::exists("dcosth")) dcosth = Config::get<double>("dcosth");
  if (Config::exists("nth")) {
    int nth = Config::get<int>("nth");
    dcosth = 2./nth;
  }
  double sum(0);
  for (double costh(-1. + dcosth / 2.); costh < 1.; costh += dcosth) {
    sum += dsig_dcosth(costh);
  }
  return sum * dcosth;
}
