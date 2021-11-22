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

#include "BornTerms.hpp"
#include "Isospin.hpp"
#include "Vrancx.hpp"
#include "oldModel.hpp"
#include "wavefunc.hpp"

double formfactorRNpi(string resonance, double m2) {
  if (Config::exists("noRNpiFF")) return 1;
  double m = (m2 > 0) ? sqrt(m2) : 0;
  double mN = Config::get<double>("Nucleon.mass");
  double mpi = Config::get<double>("pi_pm.mass");
  double mR = Config::get<double>(resonance + ".mass");
  double Gamma = Config::get<double>(resonance + ".width");
  int l = Config::get<double>(resonance + ".l");
  double delta2 = pow(mR - mN - mpi, 2) + Gamma * Gamma / 4.;
  double q0 = momentum(mR, mN, mpi);
  double q = (m > mN + mpi) ? momentum(m, mN, mpi) : 0;
  // return sqrt(mR / m) *
  //       pow((q0 * q0 + delta2) / (q * q + delta2), (l + 1.) / 2.);
  // double ret= sqrt(mR / m) *
  //       pow((q0 * q0 + delta2) / (q * q + delta2), (l + 1.) / 2.);
  double ret = pow((q0 * q0 + delta2) / (q * q + delta2),
                   (l + 1.) / 2.);  // correspond to Deniz
  if (isinf(ret)) {
    cerr << "inf in formfactorRNpi" << endl;
    PR(sqrt(mR / m));
    PR(pow((q0 * q0 + delta2) / (q * q + delta2), (l + 1.) / 2.));
    exit(0);
  }
  // PR(ret);
  return ret;
}

double uchCutoff(string resonance, double q) {
  double Lambda_u = 0.3*GeV;
  if (isSet("Lambda_u")) {
    Lambda_u = Config::get<double>("Lambda_u");
  }
  halfint J=Config::get<halfint>(resonance+".spin");
  return pow(Lambda_u/(Lambda_u+q*q),1.*J);
}

DiracMatrix vertexRNpi(string resonance, FourVector pR, FourVector pN,
                       FourVector q, uint muR1, uint muR2) {
  halfint spin = Config::get<halfint>(resonance + ".spin");
  int parity = Config::get<halfint>(resonance + ".parity");
  double mR = Config::get<double>(resonance + ".mass");
  double g = Config::get<double>(resonance + (isSet("oldModel") ? ".g0_old" : ".g0"));
  double FF = formfactorRNpi(resonance, pR * pR);
  double isofac = sqrt(2);
  if (spin == half) {
    return isofac * FF * vertex1hNpi(g, spin * parity, q);
  } else if (spin == 3 * half) {
    if (isSet("oldModel")) {
      return isofac * FF * vertex3hNpi_old(g, spin * parity, muR1, pR, q);
    }
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
  double g = Config::get<double>(resonance + (isSet("oldModel") ? ".g0_old" : ".g0"));
  double FF = formfactorRNpi(resonance, pR * pR);

  double isofac(0);
  if (isospin == half) {
    isofac = isospin_1h1h1(QR, -QN, -Qpi);
  } else if (isospin == 3 * half) {
    if (pR.future()) {  // incoming resonance
      isofac = isospin_3h1h1(QR, -QN, -Qpi);
    } else {  // outgoing resonance
      isofac = isospin_3h1h1(-QR, QN, Qpi);
    }
  } else {
    cerr << "vertexRNpi: isospin= " << isospin << " not implemented" << endl;
  }

  if (spin == half) {
    return isofac * FF * vertex1hNpi(g, spin * parity, q);
  } else if (spin == 3 * half) {
    if (isSet("oldModel")) {
      return isofac * FF * vertex3hNpi_old(g, spin * parity, muR1, pR, q);
    }
    return isofac * FF * vertex3hNpi(g, spin * parity, muR1, pR, q);
  } else if (spin == 5 * half) {
    return isofac * FF * vertex5hNpi(g, spin * parity, muR1, muR2, pR, q);
  } else {
    cerr << "vertexRNpi: spin-parity " << spin << ((parity > 0) ? "+" : "-")
         << " not implemented" << endl;
    exit(0);
  }
}

DiracMatrix vertexRNgamma(string resonance, FourVector pR, int QR,
                          FourVector pN, FourVector k, uint mu, uint muR1,
                          uint muR2) {
  halfint spin = Config::get<halfint>(resonance + ".spin");
  int parity = Config::get<halfint>(resonance + ".parity");
  double mR = Config::get<double>(resonance + ".mass");
  double g(0);
  if (QR == 0) {
    g = Config::get<double>(resonance + (isSet("oldModel") ? ".gngamma_old" : ".gngamma"));
  } else {
    g = Config::get<double>(resonance + (isSet("oldModel") ? ".gpgamma_old" : ".gpgamma"));
  }
  if (spin == half) {
    return vertex1hNgamma(g, spin * parity, pR, mu, k);
  } else if (spin == 3 * half) {
    if (isSet("oldModel")) {
      return vertex3hNgamma_old(g, spin * parity, muR1, pR, mu, k);
    }
    return vertex3hNgamma(g, 0, 0, spin * parity, muR1, pR, mu, k);
  } else if (spin == 5 * half) {
    return vertex5hNgamma(g, 0, 0, spin * parity, muR1, muR2, pR, mu, k);
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
  if (m < mN + mpi) return 0;
  double mR = Config::get<double>(resonance + ".mass");
  int l = Config::get<double>(resonance + ".l");
  double q0 = momentum(mR, mN, mpi);
  double q = momentum(m, mN, mpi);
  double FF = formfactorRNpi(resonance, m * m);
  return Gamma0 * pow(q / q0, 2. * l + 1.) * FF * FF;
}

dcomplex BreitWigner(string resonance, double s) {
  if (Config::exists("noBW")) return 1;
  double srt = (s > 0) ? sqrt(s) : 0;
  double mR(Config::get<double>(resonance + ".mass"));
  double Gamma = resonanceWidth(resonance, srt);
  return 1. / (s - mR * mR + i_ * srt * Gamma);
}

DiracMatrix propR(string resonance, FourVector p, uint muR1, uint nuR1,
                  uint muR2, uint nuR2) {
  halfint spin = Config::get<halfint>(resonance + ".spin");
  double mR = Config::get<double>(resonance + ".mass");
  double s = p * p;
  double srt = (s > 0) ? sqrt(s) : 0;
  double Gamma = resonanceWidth(resonance, srt);
  // PR(BreitWigner(resonance,s));
  if (spin == half) {
    return i_ * pro1half(p, mR) * BreitWigner(resonance, s);
  } else if (spin == 3 * half) {
    return i_ * pro1half(p, mR) * P3h(p, mR, muR1, nuR1) *
           BreitWigner(resonance, s);
  } else if (spin == 5 * half) {
    return i_ * pro1half(p, mR) * P5h(p, mR, muR1, nuR1, muR2, nuR2) *
           BreitWigner(resonance, s);
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
  int npol = 2 * JR + 1;
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
  return 1. / (8. * pi_) * kin.pabs() / (M * M) * 1. / npol * MSQR;
}

double widthRNpi(string resonance) {
  double M = Config::get<double>(resonance + ".mass");
  return widthRNpi(resonance, M);
}

double widthRNpi(string resonance, int QR, int QN, int Qpi, double M) {
  assert(QR == QN + Qpi);
  // calculate kinematics:
  halfint JR = Config::get<halfint>(resonance + ".spin");
  int npol = 2 * JR + 1;
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
    } else if (JR == 5 * half) {
      u_5h uR(pR);
      for (halfint laR(5 * half); laR >= -5 * half; laR--) {
        dcomplex helamp(0);
        for (uint mu1(0); mu1 < 4; mu1++) {
          for (uint mu2(0); mu2 < 4; mu2++) {
            helamp +=
                ubarN(0, laN) *
                vertexRNpi(resonance, pR, QR, -pN, -QN, -ppi, -Qpi, mu1, mu2) *
                uR(mu1, mu2, laR) * sign_(mu1) * sign_(mu2);
          }
        }
        MSQR += real(helamp * conj(helamp));
      }
    } else {
      cerr << "widthRNpi: JR = " << JR << " not implemented" << endl;
      exit(0);
    }
  }
  return 1. / (8. * pi_) * kin.pabs() / (M * M) * 1. / npol * MSQR;
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

double widthRNgamma(string resonance, int QR, double M) {
  halfint isospin = Config::get<halfint>(resonance + ".isospin");
  assert(QR <= isospin + half and -isospin + half <= QR);
  if (QR > 1 or QR < 0) return 0;
  // calculate kinematics:
  halfint JR = Config::get<halfint>(resonance + ".spin");
  int npol = 2 * JR + 1;
  double mN = Config::get<double>("Nucleon.mass");

  Kinema2 kin(M, mN, 0);
  FourVector pR = kin.P();
  FourVector pN = kin.p1();
  FourVector k = kin.p2();
  ubar_ ubarN(half, pN);

  double MSQR(0);
  for (halfint laN : {half, -half}) {
    if (JR == half) {
      u_1h uR(pR);
      for (halfint laR : {half, -half}) {
        for (uint mu(0); mu < 4; mu++) {
          dcomplex helamp_mu = ubarN(0, laN) *
                               vertexRNgamma(resonance, pR, QR, pN, k, mu) *
                               uR(laR);
          MSQR += -real(helamp_mu * conj(helamp_mu)) * sign_(mu);
        }
      }
    } else if (JR == 3 * half) {
      u_3h uR(pR);
      for (halfint laR : {3 * half, half, -half, -3 * half}) {
        for (uint mu(0); mu < 4; mu++) {
          dcomplex helamp_mu(0);
          for (uint nu(0); nu < 4; nu++) {
            helamp_mu += ubarN(0, laN) *
                         vertexRNgamma(resonance, pR, QR, pN, k, mu, nu) *
                         uR(nu, laR) * sign_(nu);
          }
          MSQR += -real(helamp_mu * conj(helamp_mu)) * sign_(mu);
        }
      }
    } else if (JR == 5 * half) {
      u_5h uR(pR);
      for (halfint laR(5 * half); laR >= -5 * half; laR--) {
        for (uint mu(0); mu < 4; mu++) {
          dcomplex helamp_mu(0);
          for (uint nu1(0); nu1 < 4; nu1++) {
            for (uint nu2(0); nu2 < 4; nu2++) {
              helamp_mu +=
                  ubarN(0, laN) *
                  vertexRNgamma(resonance, pR, QR, pN, k, mu, nu1, nu2) *
                  uR(nu1, nu2, laR) * sign_(nu1) * sign_(nu2);
            }
          }
          MSQR += -real(helamp_mu * conj(helamp_mu)) * sign_(mu);
        }
      }
    } else {
      cerr << "widthRNpi: JR = " << JR << " not implemented" << endl;
      exit(0);
    }
  }
  return 1. / (8. * pi_) * kin.pabs() / (M * M) * 1. / npol * MSQR;
}

double widthRNgamma(string resonance, int QR) {
  double M = Config::get<double>(resonance + ".mass");
  return widthRNgamma(resonance, QR, M);
}

// old calculation of Born terms (from 2012 paper) - isospin not implemented:
DiracMatrix BornTerms_old(FourVector pi, FourVector pf, FourVector q, uint mu) {
  double F1(1.);
  double F2(1.);
  double F3(1.);
  double Fhat(1.);
  MultiArray<DiracMatrix> M_(idx(_1, _5), idx_lor);
  DiracMatrix Tgamma;
  FourVector k = pi + q - pf;
  double s = (pi + q) * (pi + q);
  double u = (pi - k) * (pi - k);
  double t = (pi - pf) * (pi - pf);
  static const double f_NNpi = Config::get<double>("f_NNpi");
  static const double mpi = Config::get<double>("pi_pm.mass");
  static const double mn = Config::get<double>("Nucleon.mass");
  static const double mrho = Config::get<double>("rho.mass");
  static const double kan = Config::get<double>("ka_nngamma");
  static const double kap = Config::get<double>("ka_ppgamma");
  static const double kar = Config::get<double>("ka_NNrho");
  static const double alpha(1 / 137.036);
  static const double e(sqrt(4. * pi_ * alpha));
  static const double fac = e * f_NNpi / mpi * sqrt(2);
  if (Config::exists("cutoffBorn")) {
    double LBorn = 0.63;  // from Zetenyi, Wolf PRC 2012
    if (Config::exists("LBorn")) {
      LBorn = Config::get<double>("LBorn");
    }

    F1 = 1. / (1. + POW<2>(s - mn * mn) / POW<4>(LBorn));
    F2 = 1. / (1. + POW<2>(u - mn * mn) / POW<4>(LBorn));
    F3 = 1. / (1. + POW<2>(t - mpi * mpi) / POW<4>(LBorn));
    Fhat = F1 + F2 + F3 - F1 * F2 - F1 * F3 - F2 * F3 + F1 * F2 * F3;
  }
  for (halfint mu(0); mu < 4; mu++) {
    M_(_1, mu) = gamma5_ * (gamma_(mu) * gamma_(k) - k(mu) * gamma_unit);
    M_(_2, mu) = gamma5_ / 2. *
                 ((2. * pi(mu) - k(mu)) * (2. * (q * k) - k * k) -
                  (2. * q(mu) - k(mu)) * (2. * (pi * k) - k * k));
    M_(_3, mu) = gamma5_ / 2. *
                 (gamma_(mu) * (2. * (pf * k) + k * k) -
                  (2. * pf(mu) + k(mu)) * gamma_(k));
    M_(_4, mu) = gamma5_ / 2. *
                 (gamma_(mu) * (2. * (pi * k) - k * k) -
                  (2. * pi(mu) - k(mu)) * gamma_(k));
  }
  dcomplex Agamma[5];
  Agamma[1] = -fac * (1. / (2. * mn) * (F1 * kan + F2 * kap) +
                      2. * mn * F2 / (u - mn * mn) * (1. + kap) +
                      2. * mn * F1 / (s - mn * mn) * kan);
  Agamma[2] = fac * 4. * mn * Fhat / (t - mpi * mpi) / (u - mn * mn);
  Agamma[3] = fac * 2. * kan * F1 / (s - mn * mn);
  Agamma[4] = fac * 2. * kap * F2 / (u - mn * mn);
  Tgamma = gamma_null;
  for (int j(1); j <= 5; j++) {
    Tgamma += Agamma[j] * M_(_1 * j, mu);
  }
  return Tgamma;
}

pionPhotoprodTest::pionPhotoprodTest(double srt, int Qpi, int Qf)
    : srt(srt),
      Qpi(Qpi),
      Qf(Qf),
      Qi(Qf + Qpi),
      mR(Config::get<double>("N1440.mass")),
      mN(Config::get<double>("Nucleon.mass")),
      mp(Config::get<double>("proton.mass")),
      mn(Config::get<double>("neutron.mass")),
      mNi((Qi==0) ? mn : mp),
      mNf((Qf==0) ? mn : mp),
      mpipm(Config::get<double>("pi_pm.mass")),
      mpi0(Config::get<double>("pi_0.mass")),
      mpi((Qpi==0) ? mpi0 : mpipm),
      KINin(srt, 0, mNi),
      KINout(srt, mpi, mNf) {}

double pionPhotoprodTest::MSQR_numeric(double costh) {
  bool s_channel(false);
  bool u_channel(false);
  if (isSet("sch")) s_channel = true;
  if (isSet("uch")) u_channel = true;
  if (not(s_channel or u_channel)) s_channel = u_channel = true;
  int Qs = Qi;
  int Qu = Qf;
  FourVector pi = KINin.p2(1);
  FourVector k = KINin.p1(1);
  FourVector pf = KINout.p2(costh);
  FourVector q = KINout.p1(costh);
  FourVector ps = pi + k;
  FourVector pu = pf - k;
  BornTerms BT(pi, Qs, -q, -Qpi, k);
  /*
  if (costh>0.99) {
    for (int i(1); i<=5; i++) {
      cerr << "A" << i << " = " << BT.Agamma[i] << endl;
    }
  }
  */
  MultiArray<DiracMatrix> T(idx_lor);
  for (uint mu(0); mu < 4; mu++) {
    T(mu) = gamma_null;
    if (isSet("Born")) T(mu) += BT.Tgamma_all(mu);
    if (isSet("oldBorn")) {
      T(mu) += BornTerms_old(pi, pf, -q, mu);
    }
    for (string resonance : {"N1440", "N1535", "N1650", "R1hp", "R1hm"}) {
      if (isSet(resonance)) {
        if (s_channel) {
          T(mu) += vertexRNpi(resonance, ps, Qs, -pf, -Qf, -q, -Qpi) *
                   propR(resonance, ps) *
                   vertexRNgamma(resonance, -ps, -Qs, pi, k, mu);
        }
        if (u_channel) {
          double uchCut = uchCutoff(resonance,KINout.pabs());
          T(mu) += vertexRNgamma(resonance, pu, Qu, -pf, k, mu) *
                   propR(resonance, pu) *
                   vertexRNpi(resonance, -pu, -Qu, pi, Qi, -q, -Qpi) * uchCut;
          /*
          cerr << "mu = " << mu << "   --------------------------" << endl;
          PR(pu); PR(pu*pu);
          PR(vertexRNgamma(resonance, pu, Qu, -pf, k, mu));
          PR(propR(resonance, pu));
          PR(vertexRNpi(resonance, -pu, -Qu, pi, Qi, -q, -Qpi));
          PR(T(mu));
          //*/
        }
      }
    }
    for (uint nu1(0); nu1 < 4; nu1++) {
      for (uint nu2(0); nu2 < 4; nu2++) {
        for (string resonance : {"D1232", "N1520", "R3hp", "R3hm"}) {
          if (isSet(resonance)) {
            if (s_channel) {
              T(mu) += vertexRNpi(resonance, ps, Qs, -pf, -Qf, -q, -Qpi, nu1) *
                       sign_(nu1) * propR(resonance, ps, nu1, nu2) *
                       sign_(nu2) *
                       vertexRNgamma(resonance, -ps, -Qs, pi, k, mu, nu2);
            }
            if (u_channel) {
              double uchCut = uchCutoff(resonance,KINout.pabs());
              T(mu) += vertexRNgamma(resonance, pu, Qu, -pf, k, nu1) *
                       sign_(nu1) * propR(resonance, pu, nu1, nu2) *
                       sign_(nu2) *
                       vertexRNpi(resonance, -pu, -Qu, pi, Qi, -q, -Qpi, nu2) * uchCut;
            }
          }
        }
      }
    }
    for (uint nu1(0); nu1 < 4; nu1++) {
      for (uint nu2(0); nu2 < 4; nu2++) {
        for (uint nu1p(0); nu1p < 4; nu1p++) {
          for (uint nu2p(0); nu2p < 4; nu2p++) {
            for (string resonance : {"N1675", "N1780", "R5hp", "R5hm"}) {
              if (isSet(resonance)) {
                if (s_channel) {
                  T(mu) +=
                      vertexRNpi(resonance, ps, Qs, -pf, -Qf, -q, -Qpi, nu1,
                                 nu2) *
                      sign_(nu1) * sign_(nu2) *
                      propR(resonance, ps, nu1, nu2, nu1p, nu2p) * sign_(nu1p) *
                      sign_(nu2p) *
                      vertexRNgamma(resonance, -ps, -Qs, pi, k, mu, nu1p, nu2p);
                }
                if (u_channel) {
                  double uchCut = uchCutoff(resonance,KINout.pabs());
                  T(mu) += vertexRNgamma(resonance, pu, Qu, -pf, k, nu1, nu2) *
                           sign_(nu1) * sign_(nu2) *
                           propR(resonance, pu, nu1, nu2, nu1p, nu2p) *
                           sign_(nu1p) * sign_(nu2p) *
                           vertexRNpi(resonance, -pu, -Qu, pi, Qi, -q, -Qpi,
                                      nu1p, nu2p) * uchCut;
                }
              }
            }
          }
        }
      }
    }
  }
  dcomplex MSQR(0);
  DiracMatrix GG = gamma_null;
  for (uint mu(0); mu < 4; mu++) {
    GG += -proN(pf) * T(mu) * proN(pi) * adj(T(mu)) * sign_(mu);
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
  // double pin_abs = KINin.p1().spacial().abs();
  // double pout_abs = KINout.p1().spacial().abs();
  double pin_abs = KINin.pabs();
  double pout_abs = KINout.pabs();
  // cerr << "phase space: " << setw(15) << mub(1. / (32. * pi_ * srt * srt) *
  // pout_abs / pin_abs * 1. / npol) << endl; cerr << "MSQR:        " <<
  // setw(15)
  // << MSQR_numeric(costh) << endl;
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
  /*
  cout << "decay width of D(1232):" << endl;
  cout << "D++ : " << widthRNpi("D1232", 2) << endl;
  cout << "D+  : " << widthRNpi("D1232", 1) << endl;
  cout << "D0  : " << widthRNpi("D1232", 0) << endl;
  cout << "D-  : " << widthRNpi("D1232", 2) << endl;
  cout << "decay width of N(1520):" << endl;
  cout << "N+ : " << widthRNpi("N1520", 1) << endl;
  cout << "N0 : " << widthRNpi("N1520", 0) << endl;
  */
  //*
  if (isSet("fixCouplings")) {
    cout << "fixing coupling constants:" << endl;
    for (string resonance : {"D1232", "N1520", "N1440", "N1535", "N1650", "N1675",
                            "N1680", "N1700", "N1710", "D1600", "D1620"}) {
      double g0 = Config::get<double>(resonance + ".g0");
      double gngamma = Config::get<double>(resonance + ".gngamma");
      double gpgamma = Config::get<double>(resonance + ".gpgamma");
      double Gamma = Config::get<double>(resonance + ".width");
      double BRNpi = Config::get<double>(resonance + ".BNpi");
      double BRngamma = Config::get<double>(resonance + ".Bngamma");
      double BRpgamma = Config::get<double>(resonance + ".Bpgamma");
      double GNpi = Gamma * BRNpi;
      double Gngamma = Gamma * BRngamma;
      double Gpgamma = Gamma * BRpgamma;
      double GNpi_calc = widthRNpi(resonance, 1);
      double Gngamma_calc = widthRNgamma(resonance, 0);
      double Gpgamma_calc = widthRNgamma(resonance, 1);
      double fac_Npi = sqrt(GNpi / GNpi_calc);
      double fac_ngamma = sqrt(Gngamma / Gngamma_calc);
      double fac_pgamma = sqrt(Gpgamma / Gpgamma_calc);
      cout << resonance << " -> N+pi   - g0: " << setw(12) << g0 << " -> "
          << setw(12) << g0 * fac_Npi << " ( width: " << setw(12) << GNpi_calc
          << " -> " << setw(12) << GNpi << ")" << endl;
      cout << resonance << " -> n+gamma - g: " << setw(12) << gngamma << " -> "
          << setw(12) << gngamma * fac_ngamma << " ( width: " << setw(12)
          << Gngamma_calc << " -> " << setw(12) << Gngamma << ")" << endl;
      cout << resonance << " -> p+gamma - g: " << setw(12) << gpgamma << " -> "
          << setw(12) << gpgamma * fac_pgamma << " ( width: " << setw(12)
          << Gpgamma_calc << " -> " << setw(12) << Gpgamma << ")" << endl;
    }
  }
  //*/
  /*
   cout << "# decay widths" << endl;
   cout << "#" << setw(8) << "JR" << setw(15) << "width" << endl;
   for (string resonance : {"R1hp", "R1hm", "R3hp", "R3hm", "R5hp", "R5hm"}) {
     cout << setw(8) << resonance << setw(15) << widthRNpi(resonance, 1)
          << setw(15) << sqrt(0.2 / widthRNpi(resonance, 1)) << endl;
   }

   cout << "D(1232) width:" << widthRNpi("D1232") << endl;
   */
  double srt = 1.5;
  double costh = 0.904827;
  pionPhotoprodTest PP1n(srt, 1, 0);
  double msqr = PP1n.MSQR_numeric(costh);
  // PR(msqr);
  /*
  for (double costh(-1.); costh < 1.; costh += 0.1) {
    cout << setw(9) << costh << setw(15) << mub(PP1n.diffsig_numeric(costh)) <<
  endl;
  }
  */
  cout << "#" << setw(9) << "sqrt(s)" << setw(15) << "elab" << setw(15)
       << "sigtot pi0p" << setw(15) << "sigtot pi-p" << setw(15)
       << "sigtot pi+n";
  /*
  cout << setw(15) << "R1hp" << setw(30) << "R1hm" << setw(30) << "R3hp"
       << setw(30) << "R3hm";
  */
  // << setw(30) << "Breit-Wigner"
  cout << endl;
  // cout << '#' << setw(150) << "width" << endl;
  double dsrt(0.1 * GeV);
  if (Config::exists("dsrt")) dsrt = Config::get<double>("dsrt");
  for (double srt(1.1); srt < 2.1; srt += dsrt) {
    double mN = Config::get<double>("Nucleon.mass");
    double elab = (srt * srt - mN * mN) / (2. * mN);
    //*
    pionPhotoprodTest PPpi0p(srt, 0, 1);
    pionPhotoprodTest PPpimp(srt, -1, 1);
    pionPhotoprodTest PPpipn(srt, 1, 0);
    auto ppi0p = mub(PPpi0p.sigtot_numeric());
    auto ppimp = mub(PPpimp.sigtot_numeric());
    auto ppipn = mub(PPpipn.sigtot_numeric());

    cout << setw(10) << srt << setw(15) << elab << setw(15) << ppi0p << setw(15)
         << ppimp << setw(15) << ppipn;
    /*
    cout << setw(15) << resonanceWidth("R1hp", srt) << setw(15)
         << widthRNpi("R1hp", 1, srt) << setw(15) << resonanceWidth("R1hm", srt)
         << setw(15) << widthRNpi("R1hm", 1, srt) << setw(15)
         << resonanceWidth("R3hp", srt) << setw(15) << widthRNpi("R3hp", 1, srt)
         << setw(15) << resonanceWidth("R3hm", srt) << setw(15)
         << widthRNpi("R3hm", 1, srt);
    //  << BreitWigner("N1520", srt*srt) << setw(15)
    //  << POW<2>(abs(BreitWigner("N1520", srt*srt)))
    */
    cout << endl;
    //*/
  }
  return 1;
}