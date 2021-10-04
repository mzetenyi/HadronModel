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

DiracMatrix vertexRNpi(FourVector pR, FourVector pN, FourVector q) {
  double isofac = sqrt(2);
  double mpi = Config::get<double>("pi_pm.mass");
  double g = Config::get<double>("N1440.g0");
  return isofac*g/mpi * gamma5_ * gamma_(q);
}

DiracMatrix vertexRNgamma(FourVector pR, FourVector pN, FourVector k, uint mu) {
  double mrho = Config::get<double>("rho.mass");
  double g = Config::get<double>("N1440.gngamma");
  return i_*g/(2.*mrho) * (gamma_(k) * gamma_(mu) - gamma_(mu) * gamma_(k));
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

double N1440width(double srt){
  if (Config::exists("noRwidth")) return 0;
  if (Config::exists("constRwidth")) {
    return Config::get<double>("N1440.width");
  }
  //cerr << "N1440 width not implemented, using 0" << endl;
  return 0;
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

dcomplex BreitWigner(FourVector p, double m, double Gamma) {
  if (Config::exists("noBW")) return 1;
  return 1./(p*p - m*m + i_*sqrt(p*p)*Gamma);
}

DiracMatrix propR(FourVector p) {
  double mR(Config::get<double>("N1440.mass"));
  double srt = sqrt(p*p);
  double Gamma = N1440width(srt);
  return i_*pro1half(p,mR) * BreitWigner(p,mR,Gamma);
}

DiracMatrix propR(string resonance, FourVector p, uint muR1, uint nuR1, uint muR2, uint nuR2) {
  halfint spin = Config::get<halfint>(resonance + ".spin");
  double mR = Config::get<double>(resonance + ".mass");
  double Gamma = resonanceWidth(resonance,sqrt(p*p));
  if (spin == half) {
    return i_*pro1half(p,mR) * BreitWigner(p,mR,Gamma);
  } else if (spin == 3*half) {
    return i_*P3h(p,mR,muR1,nuR1) * BreitWigner(p,mR,Gamma);
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
        /*
  PR(mR);
  PR(mN);
  PR(mpi);
  */
}

double pionPhotoprodTest::MSQRraw_numeric(double costh) {
  FourVector pi = KINin.p1(1);
  FourVector k = KINin.p2(1);
  FourVector pf = KINout.p1(costh);
  FourVector q = KINout.p2(costh);
  FourVector p = pi + k;

/*
  PR(pi);
  PR(k);
  PR(pf);
  PR(q);

  PR(mpi * mpi);
  PR(q * q);
  PR(mN * mN);
  PR(pi * pi);
  PR(srt * srt);
  PR(p * p);
  PR(k * k);
  PR(pf * pf);
*/

  MultiArray<DiracMatrix> T(idx_lor);
  for (uint mu(0); mu < 4; mu++) {
    if (isSet("N1440")) T(mu) += vertexRNpi("N1440",p, -pf, -q) * propR(p) * vertexRNgamma(p, pi, k, mu);
    for (uint mu(0); mu < 4; mu++) {
      if (isSet("N1520")) T(mu) += vertexRNpi("N1520",p, -pf, -q) * propR(p) * vertexRNgamma(p, pi, k, mu);

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

double pionPhotoprodTest::MSQRraw_analytic(double costh) {
  FourVector pi = KINin.p1(1);
  FourVector k = KINin.p2(1);
  FourVector pf = KINout.p1(costh);
  FourVector q = KINout.p2(costh);
  FourVector p = pi + k;

/*
  PR(pi);
  PR(k);
  PR(pf);
  PR(q);
*/

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

double pionPhotoprodTest::diffsig_analytic(double costh) {
  int npol(4);
  double pin_abs = KINin.p1().spacial().abs();
  double pout_abs = KINout.p1().spacial().abs();
  return 1./(32.*pi_*srt*srt) * pout_abs/pin_abs * 1./npol * MSQRraw_analytic(costh);
}

double pionPhotoprodTest::diffsig_numeric(double costh) {
  int npol(4);
  double pin_abs = KINin.p1().spacial().abs();
  double pout_abs = KINout.p1().spacial().abs();
  return  1./(32.*pi_*srt*srt) * pout_abs/pin_abs * 1./npol * MSQRraw_numeric(costh);
}

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
  double srt = 1.5 * GeV;
  if (Config::exists("srt")) {
    srt = Config::get<double>("srt");
  }
  double costh = 0.5;
  if (Config::exists("costh")) {
    costh = Config::get<double>("costh");
  }
  pionPhotoprodTest PPT(srt);
  double MSQR_numeric = PPT.MSQRraw_numeric(costh);
  double MSQR_analytic = PPT.MSQRraw_analytic(costh);
  PR(srt);
  PR(costh);
  PR(MSQR_analytic);
  PR(MSQR_numeric);
  double sigtot_numeric = PPT.sigtot_numeric();
  double sigtot_analytic = PPT.sigtot_analytic();
  PR(sigtot_numeric);
  PR(mub(sigtot_numeric));
  //PR(sigtot_analytic)

  for (double srt(1.2); srt<2.1; srt+=0.1) {
    pionPhotoprodTest PPT(srt);
    cout << setw(10) << srt << setw(15) << mub(PPT.sigtot_numeric()) << endl;
  }
  
  return 1;
}