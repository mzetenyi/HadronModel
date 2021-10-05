#include "pionPhotoprodTest.hpp"

#include "Config.hpp"
#include "units.hpp"
using namespace units_GeV;

#include "Spinors.hpp"
#include "utils.hpp"
using namespace Spinors;

#include "Vectors.hpp"
using namespace Vectors;

#include "MultiArray.hpp"
#include <iostream>
#include <iomanip>
using namespace std;

#include "Vrancx.hpp"

double formfactorRNpi(string resonance, double m) {
  if (Config::exists("noRNpiFF")) return 1;
  double mN = Config::get<double>("Nucleon.mass");
  double mpi = Config::get<double>("pi_pm.mass");
  double mR = Config::get<double>(resonance+".mass");
  double Gamma = Config::get<double>(resonance+".width");
  int l = Config::get<double>(resonance+".l");
  double delta2 = pow(mR-mN-mpi,2) + Gamma*Gamma/4.;
  double q0 = momentum(mR,mN,mpi);
  double q = momentum(m,mN,mpi);
  return sqrt(mR/m) * pow((q0*q0+delta2)/(q*q+delta2),(l+1.)/2.);
}

DiracMatrix vertexRNpi(string resonance, FourVector pR, FourVector pN, FourVector q, uint muR1, uint muR2) {
  halfint spin = Config::get<halfint>(resonance + ".spin");
  int parity = Config::get<halfint>(resonance + ".parity");
  double mR = Config::get<double>(resonance + ".mass");
  double g = Config::get<double>(resonance + ".g0");
  double FF = formfactorRNpi(resonance,sqrt(pR*pR));
  double isofac = sqrt(2);
  if (spin==half) {
    return isofac * FF * vertex1hNpi(g,spin*parity,q);
  } else if (spin == 3*half) {
    return isofac * FF * vertex3hNpi(g,spin*parity,muR1,pR,q);
  } else {
    cerr << "vertexRNpi: spin-parity " << spin << ((parity>0) ? "+" : "-") << " not implemented" << endl;
    exit(0);
  }
}

DiracMatrix vertexRNgamma(string resonance, FourVector pR, FourVector pN, FourVector k, uint mu, uint muR1, uint muR2) {
  halfint spin = Config::get<halfint>(resonance + ".spin");
  int parity = Config::get<halfint>(resonance + ".parity");
  double mR = Config::get<double>(resonance + ".mass");
  double g = Config::get<double>(resonance + ".gngamma");
  if (spin == half) {
    return vertex1hNgamma(g,spin*parity,pR,mu,k);
  } else if (spin == 3*half) {
    return vertex3hNgamma(g,0,0,spin*parity,muR1,pR,mu,k);
  } else {
    cerr << "vertexRNgamma: spin-parity " << spin << ((parity>0) ? "+" : "-") << " not implemented" << endl;
    exit(0);
  }
}

DiracMatrix pro1half(FourVector p, double m) {
  return gamma_(p) + m * gamma_unit;
}

double resonanceWidth(string resonance, double m) {
  if (Config::exists("noRwidth")) return 0;
  double Gamma0 = Config::get<double>(resonance + ".width");
  if (Config::exists("constRwidth")) {
    return Gamma0;
  }
  double mN = Config::get<double>("Nucleon.mass");
  double mpi = Config::get<double>("pi_pm.mass");
  double mR = Config::get<double>(resonance+".mass");
  int l = Config::get<double>(resonance+".l");
  double q0 = momentum(mR,mN,mpi);
  double q = momentum(m,mN,mpi);
  double FF = formfactorRNpi(resonance,m);
  return Gamma0 * pow(q/q0,2.*l+1.) * FF*FF;
}

dcomplex BreitWigner(string resonance, double srt) {
  if (Config::exists("noBW")) return 1;
  double mR(Config::get<double>(resonance+".mass"));
  double Gamma = resonanceWidth(resonance,srt);
  return 1./(srt*srt - mR*mR + i_*srt*Gamma);
}

DiracMatrix propR(string resonance, FourVector p, uint muR1, uint nuR1, uint muR2, uint nuR2) {
  halfint spin = Config::get<halfint>(resonance + ".spin");
  double mR = Config::get<double>(resonance + ".mass");
  double Gamma = resonanceWidth(resonance,sqrt(p*p));
  double srt = sqrt(p*p);
  if (spin == half) {
    return i_*pro1half(p,mR) * BreitWigner(resonance,srt);
  } else if (spin == 3*half) {
    return i_*pro1half(p,mR) * P3h(p,mR,muR1,nuR1) * BreitWigner(resonance,srt);
  } else {
    cerr << "propR: spin " << spin << " not implemented" << endl;
    exit(0);
  }
}


DiracMatrix proN(FourVector p) {
  double mN(Config::get<double>("Nucleon.mass"));
  return gamma_(p) + mN * gamma_unit;
}

pionPhotoprodTest::pionPhotoprodTest(double srt)
    : srt(srt),
      mR(Config::get<double>("N1440.mass")),
      mN(Config::get<double>("Nucleon.mass")),
      mpi(Config::get<double>("pi_pm.mass")),
      KINin(srt, mN, 0),
      KINout(srt, mN, mpi) {
}

double pionPhotoprodTest::MSQRraw_numeric(double costh) {
  FourVector pi = KINin.p1(1);
  FourVector k = KINin.p2(1);
  FourVector pf = KINout.p1(costh);
  FourVector q = KINout.p2(costh);
  FourVector p = pi + k;
  MultiArray<DiracMatrix> T(idx_lor);
  for (uint mu(0); mu < 4; mu++) {
    if (isSet("N1440")) T(mu) += vertexRNpi("N1440",p, -pf, -q) * propR("N1440",p) * vertexRNgamma("N1440",p, pi, k, mu);
    for (uint nu1(0); nu1 < 4; nu1++) {
      for (uint nu2(0); nu2 < 4; nu2++) {
        if (isSet("N1520")) T(mu) += vertexRNpi("N1520",p, -pf, -q, nu1) * propR("N1520",p,nu1,nu2) * 
          vertexRNgamma("N1520",-p, pi, k, mu, nu2);
      }
    }
  }
  dcomplex MSQR(0);
  DiracMatrix GG = gamma_null;
  for (uint mu(0); mu < 4; mu++) {
    GG += - proN(pf) * T(mu) * proN(pi) * adj(T(mu) * sign_(mu));
  }
  MSQR = trace(GG);
  return real(MSQR);
}
/*
double pionPhotoprodTest::MSQRraw_analytic(double costh) {
  FourVector pi = KINin.p1(1);
  FourVector k = KINin.p2(1);
  FourVector pf = KINout.p1(costh);
  FourVector q = KINout.p2(costh);
  FourVector p = pi + k;

  double pi_pf = pi * pf;
  double pf_q = pf * q;
  double pi_q = pi * q;
  double pi_k = pi * k;
  double k_p = k * p;
  double p_pf = p * pf;
  double p_q = p * q;
  double k_pf = k * pf;
  double k_q = k * q;

  double MSQR = -128 * pi_k * k_p * p_pf * POW<2>(mpi) +
                256 * pi_k * k_p * p_q * pf_q -
                128 * pi_k * k_p * mR * mN * POW<2>(mpi) +
                64 * pi_k * k_pf * POW<2>(mpi) * POW<2>(srt) -
                64 * pi_k * k_pf * POW<2>(mR) * POW<2>(mpi) -
                128 * pi_k * k_q * pf_q * POW<2>(srt) +
                128 * pi_k * k_q * pf_q * POW<2>(mR);

  return MSQR;
}
*/
/*
double pionPhotoprodTest::diffsig_analytic(double costh) {
  int npol(4);
  double pin_abs = KINin.p1().spacial().abs();
  double pout_abs = KINout.p1().spacial().abs();
  return 1./(32.*pi_*srt*srt) * pout_abs/pin_abs * 1./npol * MSQRraw_analytic(costh);
}
*/
double pionPhotoprodTest::diffsig_numeric(double costh) {
  int npol(4);
  double pin_abs = KINin.p1().spacial().abs();
  double pout_abs = KINout.p1().spacial().abs();
  return  1./(32.*pi_*srt*srt) * pout_abs/pin_abs * 1./npol * MSQRraw_numeric(costh);
}
/*
double pionPhotoprodTest::sigtot_analytic() {
  double dcosth = 0.01;
  if (Config::exists("dcosth")) dcosth = Config::get<double>("dcosth");
  if (Config::exists("nth")) {
    int nth = Config::get<int>("nth");
    dcosth = 2. / nth;
  }
  double sum(0);
  for (double costh(-1. + dcosth / 2.); costh < 1.; costh += dcosth) {
    sum += diffsig_analytic(costh);
  }
  return sum * dcosth;
}
*/
double pionPhotoprodTest::sigtot_numeric() {
  double dcosth = 0.01;
  if (Config::exists("dcosth")) dcosth = Config::get<double>("dcosth");
  if (Config::exists("nth")) {
    int nth = Config::get<int>("nth");
    dcosth = 2. / nth;
  }
  double sum(0);
  for (double costh(-1. + dcosth / 2.); costh < 1.; costh += dcosth) {
    sum += diffsig_numeric(costh);
  }
  return sum * dcosth;
}

int main(int argc, char** argv) {
  Config::load(argc, argv);
  
  cout << "#" << setw(9) << "sqrt(s)" << setw(15) << "sigtot [mub]" << setw(30) << "Breit-Wigner" << endl;
  for (double srt(1.2); srt<3.1; srt+=0.1) {
    pionPhotoprodTest PPT(srt);
    cout << setw(10) << srt << setw(15) << mub(PPT.sigtot_numeric()) 
        << setw(30) << BreitWigner("N1520",srt) << setw(15) << POW<2>(abs(BreitWigner("N1520",srt))) << endl;
  }
  
  return 1;
}