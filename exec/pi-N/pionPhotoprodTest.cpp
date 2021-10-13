#include "pionPhotoprodTest.hpp"

#include "Config.hpp"
#include "units.hpp"
using namespace units_GeV;

#include "Spinors.hpp"
#include "utils.hpp"
using namespace Spinors;

#include "Vectors.hpp"
using namespace Vectors;

#include <cassert>
#include <iomanip>
#include <iostream>

#include "MultiArray.hpp"
using namespace std;

#include "Isospin.hpp"
#include "Vrancx.hpp"
#include "wavefunc.hpp"

double formfactorRNpi(string resonance, double m) {
  if (Config::exists("noRNpiFF")) return 1;
  double mN = Config::get<double>("Nucleon.mass");
  double mpi = Config::get<double>("pi_pm.mass");
  double mR = Config::get<double>(resonance + ".mass");
  double Gamma = Config::get<double>(resonance + ".width");
  int l = Config::get<double>(resonance + ".l");
  double delta2 = pow(mR - mN - mpi, 2) + Gamma * Gamma / 4.;
  double q0 = momentum(mR, mN, mpi);
  double q = momentum(m, mN, mpi);
  return sqrt(mR / m) *
         pow((q0 * q0 + delta2) / (q * q + delta2), (l + 1.) / 2.);
}

DiracMatrix vertexRNpi(string resonance, FourVector pR, FourVector pN,
                       FourVector q, uint muR1, uint muR2) {
  halfint spin = Config::get<halfint>(resonance + ".spin");
  int parity = Config::get<halfint>(resonance + ".parity");
  double mR = Config::get<double>(resonance + ".mass");
  double g = Config::get<double>(resonance + ".g0");
  double FF = formfactorRNpi(resonance, sqrt(pR * pR));
  double isofac = sqrt(2);
  if (spin == half) {
    return isofac * FF * vertex1hNpi(g, spin * parity, q);
  } else if (spin == 3 * half) {
    return isofac * FF * vertex3hNpi(g, spin * parity, muR1, pR, q);
  } else {
    cerr << "vertexRNpi: spin-parity " << spin << ((parity > 0) ? "+" : "-")
         << " not implemented" << endl;
    exit(0);
  }
}

DiracMatrix vertexRNpi(string resonance, FourVector pR, int QR, FourVector pN,
                       int QN, FourVector q, int Qpi, uint muR1, uint muR2) {
  halfint spin = Config::get<halfint>(resonance + ".spin");
  halfint isospin = Config::get<halfint>(resonance + ".isospin");
  int parity = Config::get<halfint>(resonance + ".parity");
  double mR = Config::get<double>(resonance + ".mass");
  double g = Config::get<double>(resonance + ".g0");
  double FF = formfactorRNpi(resonance, sqrt(pR * pR));

  double isofac(0);
  if (isospin == half) {
    isofac = isospin_1h1h1(QR, -QN, -Qpi);
  } else if (isospin == 3 * half) {
    if (pR.future()) {
      isofac = isospin_3h1h1(QR, -QN, -Qpi);
    } else {
      isofac = isospin_1h3h1(QR, -QN, -Qpi);
    }
  } else {
    cerr << "vertexRNpi: isospin=" << isospin << " not implemented" << endl;
  }

  if (spin == half) {
    return isofac * FF * vertex1hNpi(g, spin * parity, q);
  } else if (spin == 3 * half) {
    return isofac * FF * vertex3hNpi(g, spin * parity, muR1, pR, q);
  } else {
    cerr << "vertexRNpi: spin-parity " << spin << ((parity > 0) ? "+" : "-")
         << " not implemented" << endl;
    exit(0);
  }
}

DiracMatrix vertexRNgamma(string resonance, FourVector pR, FourVector pN,
                          FourVector k, uint mu, uint muR1, uint muR2) {
  halfint spin = Config::get<halfint>(resonance + ".spin");
  int parity = Config::get<halfint>(resonance + ".parity");
  double mR = Config::get<double>(resonance + ".mass");
  double g = Config::get<double>(resonance + ".gngamma");
  if (spin == half) {
    return vertex1hNgamma(g, spin * parity, pR, mu, k);
  } else if (spin == 3 * half) {
    return vertex3hNgamma(g, 0, 0, spin * parity, muR1, pR, mu, k);
  } else {
    cerr << "vertexRNgamma: spin-parity " << spin << ((parity > 0) ? "+" : "-")
         << " not implemented" << endl;
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
  double mR = Config::get<double>(resonance + ".mass");
  int l = Config::get<double>(resonance + ".l");
  double q0 = momentum(mR, mN, mpi);
  double q = momentum(m, mN, mpi);
  double FF = formfactorRNpi(resonance, m);
  return Gamma0 * pow(q / q0, 2. * l + 1.) * FF * FF;
}

dcomplex BreitWigner(string resonance, double srt) {
  if (Config::exists("noBW")) return 1;
  double mR(Config::get<double>(resonance + ".mass"));
  double Gamma = resonanceWidth(resonance, srt);
  return 1. / (srt * srt - mR * mR + i_ * srt * Gamma);
}

DiracMatrix propR(string resonance, FourVector p, uint muR1, uint nuR1,
                  uint muR2, uint nuR2) {
  halfint spin = Config::get<halfint>(resonance + ".spin");
  double mR = Config::get<double>(resonance + ".mass");
  double Gamma = resonanceWidth(resonance, sqrt(p * p));
  double srt = sqrt(p * p);
  if (spin == half) {
    return i_ * pro1half(p, mR) * BreitWigner(resonance, srt);
  } else if (spin == 3 * half) {
    return i_ * pro1half(p, mR) * P3h(p, mR, muR1, nuR1) *
           BreitWigner(resonance, srt);
  } else {
    cerr << "propR: spin " << spin << " not implemented" << endl;
    exit(0);
  }
}

DiracMatrix proN(FourVector p) {
  double mN(Config::get<double>("Nucleon.mass"));
  return gamma_(p) + mN * gamma_unit;
}

double widthRNpi(string resonance, double M) {
  // calculate kinematics:
  halfint JR = Config::get<halfint>(resonance + ".spin");
  double mN = Config::get<double>("Nucleon.mass");
  double mpi(Config::get<double>("pi_pm.mass"));

  Kinema2 kin(M, mN, mpi);
  FourVector pR = kin.P();
  FourVector pN = kin.p1();
  FourVector ppi = kin.p2();
  ubar_ ubarN(half, pN);

  double MSQR(0);
  for (halfint laN : {half, -half}) {
    if (JR == half) {
      u_1h uR(pR);
      for (halfint laR : {half, -half}) {
        dcomplex helamp =
            ubarN(0, laN) * vertexRNpi(resonance, pR, -pN, -ppi) * uR(laR);
        MSQR += real(helamp * conj(helamp));
      }
    } else if (JR == 3 * half) {
      u_3h uR(pR);
      for (halfint laR : {3 * half, half, -half, -3 * half}) {
        dcomplex helamp(0);
        for (uint mu(0); mu < 4; mu++) {
          helamp += ubarN(0, laN) * vertexRNpi(resonance, pR, -pN, -ppi, mu) *
                    uR(mu, laR) * sign_(mu);
        }
        MSQR += real(helamp * conj(helamp));
      }
    } else {
      cerr << "widthRNpi: JR = " << JR << " not implemented" << endl;
      exit(0);
    }
  }
  return 1. / (8. * pi_) * kin.pabs() / (M * M) * MSQR;
}

double widthRNpi(string resonance) {
  double M = Config::get<double>(resonance + ".mass");
  return widthRNpi(resonance, M);
}

double widthRNpi(string resonance, int QR, int QN, int Qpi, double M) {
  assert(QR == QN + Qpi);
  // calculate kinematics:
  halfint JR = Config::get<halfint>(resonance + ".spin");
  double mN = Config::get<double>("Nucleon.mass");
  double mpi(Config::get<double>("pi_pm.mass"));

  Kinema2 kin(M, mN, mpi);
  FourVector pR = kin.P();
  FourVector pN = kin.p1();
  FourVector ppi = kin.p2();
  ubar_ ubarN(half, pN);

  double MSQR(0);
  for (halfint laN : {half, -half}) {
    if (JR == half) {
      u_1h uR(pR);
      for (halfint laR : {half, -half}) {
        dcomplex helamp = ubarN(0, laN) *
                          vertexRNpi(resonance, pR, QR, -pN, -QN, -ppi, -Qpi) *
                          uR(laR);
        MSQR += real(helamp * conj(helamp));
      }
    } else if (JR == 3 * half) {
      u_3h uR(pR);
      for (halfint laR : {3 * half, half, -half, -3 * half}) {
        dcomplex helamp(0);
        for (uint mu(0); mu < 4; mu++) {
          helamp += ubarN(0, laN) *
                    vertexRNpi(resonance, pR, QR, -pN, -QN, -ppi, -Qpi, mu) *
                    uR(mu, laR) * sign_(mu);
        }
        MSQR += real(helamp * conj(helamp));
      }
    } else {
      cerr << "widthRNpi: JR = " << JR << " not implemented" << endl;
      exit(0);
    }
  }
  return 1. / (8. * pi_) * kin.pabs() / (M * M) * MSQR;
}

double widthRNpi(string resonance, int QR, double M) {
  halfint isospin = Config::get<halfint>(resonance + ".isospin");
  assert(-isospin + half <= QR and QR <= isospin + half);
  double width(0);
  for (int QN : {1, 0}) {
    int Qpi = QR - QN;
    if (-1 <= Qpi and Qpi <= 1) {
      width += widthRNpi(resonance, QR, QN, Qpi, M);
    }
  }
  return width;
}

double widthRNpi(string resonance, int QR, int QN, int Qpi) {
  double M = Config::get<double>(resonance + ".mass");
  return widthRNpi(resonance, QR, QN, Qpi, M);
}

double widthRNpi(string resonance, int QR) {
  double M = Config::get<double>(resonance + ".mass");
  return widthRNpi(resonance, QR, M);
}

pionPhotoprodTest::pionPhotoprodTest(double srt, int Qpi, int QN)
    : srt(srt),
      Qpi(Qpi),
      QN(QN),
      mR(Config::get<double>("N1440.mass")),
      mN(Config::get<double>("Nucleon.mass")),
      mpi(Config::get<double>("pi_pm.mass")),
      KINin(srt, mN, 0),
      KINout(srt, mN, mpi) {}

double pionPhotoprodTest::MSQR_numeric(double costh) {
  int QR = QN + Qpi;
  FourVector pi = KINin.p1(1);
  FourVector k = KINin.p2(1);
  FourVector pf = KINout.p1(costh);
  FourVector q = KINout.p2(costh);
  FourVector p = pi + k;
  MultiArray<DiracMatrix> T(idx_lor);
  for (uint mu(0); mu < 4; mu++) {
    T(mu) = gamma_null;
    for (string resonance : {"N1440", "N1535", "N1650", "R1hp", "R1hm"}) {
      if (isSet(resonance)) {
        T(mu) += vertexRNpi(resonance, p, QR, -pf, -QN, -q, -Qpi) *
                 propR(resonance, p) * vertexRNgamma(resonance, p, pi, k, mu);
      }
    }
    for (uint nu1(0); nu1 < 4; nu1++) {
      for (uint nu2(0); nu2 < 4; nu2++) {
        for (string resonance : {"D1232", "N1520", "R3hp", "R3hm"}) {
          if (isSet(resonance)) {
            T(mu) += vertexRNpi(resonance, p, QR, -pf, -QN, -q, -Qpi, nu1) *
                     propR(resonance, p, nu1, nu2) *
                     vertexRNgamma(resonance, -p, pi, k, mu, nu2) * sign_(nu1) *
                     sign_(nu2);
          }
        }
      }
    }
  }
  dcomplex MSQR(0);
  DiracMatrix GG = gamma_null;
  for (uint mu(0); mu < 4; mu++) {
    GG += -proN(pf) * T(mu) * proN(pi) * adj(T(mu) * sign_(mu));
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
  return 1./(32.*pi_*srt*srt) * pout_abs/pin_abs * 1./npol *
MSQRraw_analytic(costh);
}
*/
double pionPhotoprodTest::diffsig_numeric(double costh) {
  int npol(4);
  double pin_abs = KINin.p1().spacial().abs();
  double pout_abs = KINout.p1().spacial().abs();
  return 1. / (32. * pi_ * srt * srt) * pout_abs / pin_abs * 1. / npol *
         MSQR_numeric(costh);
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
  cout << "decay width of D(1232):" << endl;
  cout << "D++ : " << widthRNpi("D1232", 2) << endl;
  cout << "D+  : " << widthRNpi("D1232", 1) << endl;
  cout << "D0  : " << widthRNpi("D1232", 0) << endl;
  cout << "D-  : " << widthRNpi("D1232", 2) << endl;
  cout << "decay width of N(1520):" << endl;
  cout << "N+ : " << widthRNpi("N1520", 1) << endl;
  cout << "N0 : " << widthRNpi("N1520", 0) << endl;

  cout << "fixing coupling constants:" << endl;
  for (string resonance : {"D1232", "N1520", "N1440", "N1535"}) {
    double g0 = Config::get<double>(resonance + ".g0");
    double Gamma = Config::get<double>(resonance + ".width");
    double Gamma_calc = widthRNpi(resonance, 1);
    double fac = sqrt(Gamma / Gamma_calc);
    cout << resonance << " - g0: " << g0 << " -> " << g0 * fac
         << "( width: " << Gamma_calc << " -> " << Gamma << ")" << endl;
  }
  cout << "# decay widths" << endl;
  cout << "#" << setw(8) << "JR" << setw(15) << "width" << endl;
  for (string resonance : {"R1hp", "R1hm", "R3hp", "R3hm"}) {
    cout << setw(8) << resonance << setw(15) << widthRNpi(resonance) << setw(15)
         << sqrt(0.2 / widthRNpi(resonance)) << endl;
  }

  cout << "D(1232) width:" << widthRNpi("D1232") << endl;
  cout << "#" << setw(9) << "sqrt(s)" << setw(15) << "sigtot pi0p" << setw(15)
       << "sigtot pi-p" << setw(15) << "sigtot pi+n" << setw(15) << "R1hp"
       << setw(30) << "R1hm" << setw(30) << "R3hp" << setw(30)
       << "R3hm"
       // << setw(30) << "Breit-Wigner"
       << endl;
  cout << '#' << setw(50) << "width" << endl;
  double dsrt(0.1 * GeV);
  if (Config::exists("dsrt")) dsrt = Config::get<double>("dsrt");
  for (double srt(1.1); srt < 2.1; srt += dsrt) {
    pionPhotoprodTest PPpi0p(srt, 0, 1);
    pionPhotoprodTest PPpimp(srt, -1, 1);
    pionPhotoprodTest PPpipn(srt, 1, 0);
    cout << setw(10) << srt << setw(15) << mub(PPpi0p.sigtot_numeric())
         << setw(15) << mub(PPpimp.sigtot_numeric()) << setw(15)
         << mub(PPpipn.sigtot_numeric()) << setw(15)
         << resonanceWidth("R1hp", srt) << setw(15) << widthRNpi("R1hp", srt)
         << setw(15) << resonanceWidth("R1hm", srt) << setw(15)
         << widthRNpi("R1hm", srt) << setw(15) << resonanceWidth("R3hp", srt)
         << setw(15) << widthRNpi("R3hp", srt) << setw(15)
         << resonanceWidth("R3hm", srt) << setw(15)
         << widthRNpi("R3hm", srt)
         //  << BreitWigner("N1520", srt) << setw(15)
         //  << POW<2>(abs(BreitWigner("N1520", srt)))
         << endl;
  }

  return 1;
}