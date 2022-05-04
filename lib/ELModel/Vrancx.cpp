#include "Vrancx.hpp"

#include <cassert>

#include "Config.hpp"
#include "utils.hpp"

// bool vectorTerm() {
//   static bool vterm(Config::exists("vector") or
//                     ((not Config::exists("vector")) and (not
//                     Config::exists("tensor")))
//                     );
//   //  clog << "vterm = " << vterm << endl;
//   return vterm;
// }

// bool tensorTerm() {
//   static bool tterm(Config::exists("tensor") or
//                     ((not Config::exists("vector")) and (not
//                     Config::exists("tensor")))
//                     );
//   //  clog << "tterm = " << tterm << endl;
//   return tterm;
// }

bool vectorRho() {
  static bool vRho(isSet("vector") or isSet("vectorRho") or
                   (not(isSet("tensor") or isSet("tensorRho") or
                        isSet("tensorGamma") or isSet("vectorGamma"))));
  return vRho;
}

bool tensorRho() {
  static bool tRho(isSet("tensor") or isSet("tensorRho") or
                   (not(isSet("vector") or isSet("vectorRho") or
                        isSet("vectorGamma") or isSet("tensorGamma"))));
  return tRho;
}

bool vectorGamma() {
  static bool vGamma(isSet("vector") or isSet("vectorGamma") or
                     (not(isSet("tensor") or isSet("tensorGamma") or
                          isSet("tensorRho") or isSet("vectorRho"))));
  return vGamma;
}

bool tensorGamma() {
  static bool tGamma(isSet("tensor") or isSet("tensorGamma") or
                     (not(isSet("vector") or isSet("vectorGamma") or
                          isSet("vectorRho") or isSet("tensorRho"))));
  return tGamma;
}

bool vectorOmega() {
  static bool vOmega(isSet("vector") or isSet("vectorOmega") or
                     (not(isSet("tensor") or isSet("tensorOmega"))));
  return vOmega;
}

bool tensorOmega() {
  static bool tOmega(isSet("tensor") or isSet("tensorOmega") or
                     (not(isSet("vector") or isSet("vectorOmega"))));
  return tOmega;
}

double O32(FourVector p, uint mu, uint nu, uint la) {
  return p(mu) * g_(nu, la) - p(nu) * g_(mu, la);
}

DiracMatrix O32(FourVector p, uint mu, uint la) {
  return p(mu) * gamma_(la) - g_(mu, la) * gamma_(p);
}

double O52(FourVector p, uint mu, uint nu, uint la, uint ro, uint si, uint ta) {
  return 1. / 4. *
         (+O32(p, mu, la, si) * O32(p, nu, ro, ta) +
          O32(p, mu, ro, si) * O32(p, nu, la, ta) +
          O32(p, mu, la, ta) * O32(p, nu, ro, si) +
          O32(p, mu, ro, ta) * O32(p, nu, la, si));
}

double O52(FourVector p, uint mu, uint nu, uint si, uint ta) {
  return 1. / 2. *
         (+2 * p(mu) * p(nu) * g_(si, ta) +
          (p * p) * (g_(mu, si) * g_(nu, ta) + g_(mu, ta) * g_(nu, si)) -
          p(mu) * p(si) * g_(nu, ta) - p(nu) * p(ta) * g_(mu, si) -
          p(mu) * p(ta) * g_(nu, si) - p(nu) * p(si) * g_(mu, ta));
}

double iso1hh(int Qin, int Qmes) {
  assert(Qin == 0 or Qin == 1);
  assert(Qmes == -1 or Qmes == 0 or Qmes == 1);
  if (Qmes == 0) {
    if (Qin == 0) return -1;
    return 1;
  }
  return sqrt(2);
}

double isoNNrhopi(int Qin, int Qrho, int Qpi) {
  assert(Qin == 0 or Qin == 1);
  assert(Qrho != Qpi);
  int isorhopi(1);
  if (Qrho == 1 and Qpi == -1) isorhopi = -1;
  if (Qrho == 0 and Qpi == 1) isorhopi = -1;
  if (Qrho == -1 and Qpi == 0) isorhopi = -1;
  return isorhopi * iso1hh(Qin, Qrho + Qpi);
}

dcomplex vertexrhopipi(uint mu, FourVector qp, FourVector qm) {
  static const double grho_tilde = Config::get<double>("grho_tilde");
  return i_ * grho_tilde * (qp(mu) - qm(mu));
}

dcomplex vertexgammapipi(uint mu, FourVector qp, FourVector qm) {
  return -i_ * e * (qp(mu) - qm(mu));
}

dcomplex vertexsipipi(FourVector q1, FourVector q2) {
  static const double mpi = Config::get<double>("pi_pm.mass");
  static const double gsipipi = Config::get<double>("gsipipi");
  return -i_ * gsipipi / (2 * mpi) * (q1 * q2);
}

DiracMatrix vertexNNpi(FourVector q, int Qin, int Qmes) {
  static const double mpi = Config::get<double>("pi_pm.mass");
  static const double f_NNpi = Config::get<double>("f_NNpi");
  return -f_NNpi / mpi * iso1hh(Qin, Qmes) * gamma5_ * gamma_(q);
}

DiracMatrix vertexNNsigma(double g) {
  static const double g_NNsigma = Config::get<double>("g_NNsigma");
  return i_ * g_NNsigma * gamma_unit;
}

DiracMatrix vertexNNgamma(uint mu, int Qin, FourVector k) {
  assert(Qin == 0 or Qin == 1);
  static const double mn = Config::get<double>("Nucleon.mass");
  static const double kan(Config::get<double>("ka_nngamma"));
  static const double kap(Config::get<double>("ka_ppgamma"));
  const double kappa = (Qin == 0 ? kan : kap);
  static const bool vTerm(vectorGamma());
  static const bool tTerm(tensorGamma());
  DiracMatrix ret = gamma_null;
  if (vTerm) {
    ret += -i_ * e * Qin * gamma_(mu);
  }
  if (tTerm) {
    //    ret += -i_*e * kappa/(4.*mn)*(gamma_(mu)*gamma_(k) -
    //    gamma_(k)*gamma_(mu));  // original - sign mistake c.f. Ericson-Weise
    ret += i_ * e * kappa / (4. * mn) *
           (gamma_(mu) * gamma_(k) - gamma_(k) * gamma_(mu));
  }
  return ret;
}

DiracMatrix vertexNNrho(uint mu, FourVector k, int Qin, int Qmes) {
  static const double mn = Config::get<double>("Nucleon.mass");
  static const double grho_tilde = Config::get<double>("grho_tilde");
  static const double kappa = Config::get<double>("ka_NNrho");
  static const bool vTerm(vectorRho());
  static const bool tTerm(tensorRho());
  DiracMatrix ret = gamma_null;
  if (vTerm) {
    ret += i_ * grho_tilde / 2. * iso1hh(Qin, Qmes) *
           gamma_(mu);  // only vector term
  }
  if (tTerm) {
    ret +=
        i_ * grho_tilde / 2. * iso1hh(Qin, Qmes) * kappa / (4. * mn) *
        (gamma_(mu) * gamma_(k) - gamma_(k) * gamma_(mu));  // only tensor term
  }
  return ret;
}

DiracMatrix vertexNNomega(uint mu, FourVector k) {
  static const double mn = Config::get<double>("Nucleon.mass");
  static const double g_NNomega = Config::get<double>("g_NNomega");
  static const double kappa = Config::get<double>("ka_NNomega");
  static const bool vTerm(vectorOmega());
  static const bool tTerm(tensorOmega());
  DiracMatrix ret = gamma_null;
  if (vTerm) {
    ret += i_ * g_NNomega / 2. * gamma_(mu);
  }
  if (tTerm) {
    ret += i_ * g_NNomega / 2. * kappa / (4. * mn) *
           (gamma_(mu) * gamma_(k) - gamma_(k) * gamma_(mu));
  }
  return ret;
}

DiracMatrix vertexNNrhopi(uint mu, int Qin, int Qrho, int Qpi) {
  static const double mpi = Config::get<double>("pi_pm.mass");
  static const double f_NNpi = Config::get<double>("f_NNpi");
  static const double grho_tilde = Config::get<double>("grho_tilde");
  static const double g = f_NNpi * grho_tilde;
  return g / mpi * isoNNrhopi(Qin, Qrho, Qpi) * gamma5_ * gamma_(mu);
}

DiracMatrix vertexNNgammapi(uint mu) {
  static const double mpi = Config::get<double>("pi_pm.mass");
  static const double f_NNpi = Config::get<double>("f_NNpi");
  static const double g = f_NNpi * e;
  return -g / mpi * gamma5_ * gamma_(mu);
}

DiracMatrix vertex1hNpi(double g, halfint spinParity, FourVector q) {
  const double mpi = Config::get<double>("pi_pm.mass");
  DiracMatrix ret = -g / mpi * gamma_(q);
  return ((spinParity > 0) ? gamma5_ * ret : ret);
}

DiracMatrix vertex1hNsi(double g, halfint spinParity) {
  DiracMatrix ret = -i_ * g * gamma_unit;
  return ((spinParity > 0) ? ret : ret * gamma5_);
}

DiracMatrix vertex1hNrho(double g1, halfint spinParity, FourVector pR, uint nu,
                         FourVector k) {
  const double mrho = Config::get<double>("rho.mass");
  DiracMatrix ret = g1 / mrho * sigma_(k, nu);
  return ((spinParity > 0) ? ret
                           : ((pR(0) > 0) ? -gamma5_ * ret : gamma5_ * ret));
}

DiracMatrix vertex1hNgamma(double g1, halfint spinParity, FourVector pR,
                           uint nu, FourVector k) {
  return vertex1hNrho(g1, spinParity, pR, nu, k);
}

DiracMatrix vertex3hNpi(double g, halfint spinParity, uint muR, FourVector pR,
                        FourVector q) {
  const double mpi = Config::get<double>("pi_pm.mass");
  DiracMatrix ret =
      -i_ * g / (mpi * mpi) *
      ((pR * q) * gamma_(muR) - q(muR) * gamma_(pR));  // ORIGINAL !!!!!!!!
  // DiracMatrix ret = ((pR*q)*gamma_(muR) - q(muR)*gamma_(pR));
  // cerr << "RNpi: " << g/(mpi*mpi) << endl;
  // cerr << "g, mpi: " << g << " " << mpi << endl;
  return ((spinParity > 0) ? ret : ret * gamma5_);
}

DiracMatrix vertex3hNsi(double g, halfint spinParity, uint muR, FourVector pR,
                        FourVector q) {
  const double msi = Config::get<double>("sigma.mass");
  DiracMatrix ret =
      -i_ * g / (msi * msi) *
      ((pR * q) * gamma_(muR) - q(muR) * gamma_(pR));  // ORIGINAL !!!!!!!!
  // DiracMatrix ret = g/(msi*msi) * ((pR*q)*gamma_(muR) - q(muR)*gamma_(pR));
  return ((spinParity > 0) ? ret * gamma5_ : ret);
}

DiracMatrix vertex3h3hpi(double g, halfint spinParityR, uint muR, FourVector pR,
                         halfint spinParityD, uint muD, FourVector pD,
                         FourVector q) {
  const double mpi = Config::get<double>(
      "pi_pm.mass");  // we consider only the charged pion case
  DiracMatrix Gamma = ((spinParityR == spinParityD) ? gamma5_ : gamma_unit);
  DiracMatrix ret = gamma_null;
  for (uint al(0); al < 4; al++) {
    if (pR(0) > 0) {  // incoming resonance
      ret += g / POW<3>(mpi) * O32(pD, al, muD) * Gamma * gamma_(q) *
             O32(pR, al, muR) * sign_(al);
    } else {
      ret += g / POW<3>(mpi) * O32(pR, al, muR) * Gamma * gamma_(q) *
             O32(pD, al, muD) * sign_(al);
    }
  }
  return ret;
}

DiracMatrix vertex3hNrho(double g1, double g2, double g3, halfint spinParity,
                         uint muR, FourVector pR, uint nu, FourVector k) {
  const double mrho = Config::get<double>("rho.mass");
  FourVector pN(-pR - k);
  DiracMatrix ret1(gamma_null);
  DiracMatrix ret2(gamma_null);
  DiracMatrix ret3(gamma_null);
  for (uint ro : {0, 1, 2, 3}) {
    if (g1 != 0) {
      if (pR(0) > 0) {  // incoming resonance - only valid if antibaryons are
                        // not present!
        ret1 += O32(k, ro, nu) * O32(pR, ro, muR) * sign_(ro);
      } else {  // outgiong resonance
        ret1 += O32(pR, ro, muR) * O32(k, ro, nu) * sign_(ro);
      }
    }
    if (g2 != 0) {
      ret2 += O32(pR, ro, muR) * ((pN * k) * g_(ro, nu) - k(ro) * pN(nu)) *
              sign_(ro);
    }
    if (g3 != 0) {
      ret3 +=
          O32(pR, ro, muR) * ((k * k) * g_(ro, nu) - k(ro) * k(nu)) * sign_(ro);
    }
  }
  if (spinParity > 0) {
    ret1 = ret1 * ((pR(0) > 0) ? gamma5_ : -gamma5_);
    ret2 = ret2 * gamma5_;
    ret3 = ret3 * gamma5_;
  }
  // cerr << "RNrho: " << g1/(4.*mrho*mrho) << endl;
  // cerr << "g1, mrho: " << g1 << " " << mrho << endl;
  // return ret1;
  return i_ * g1 / (4. * mrho * mrho) * ret1 -
         1. / (8. * mrho * mrho * mrho) * (g2 * ret2 + g3 * ret3);
}

DiracMatrix vertexN3hrho(double g1, double g2, double g3,
                             halfint spinParity, uint muR, FourVector pR,
                             uint nu, FourVector k) {
  const double mrho = Config::get<double>("rho.mass");
  FourVector pN(-pR - k);
  DiracMatrix ret1(gamma_null);
  DiracMatrix ret2(gamma_null);
  DiracMatrix ret3(gamma_null);
  for (uint ro : {0, 1, 2, 3}) {
    if (g1 != 0) {
      ret1 += O32(k, ro, nu) * O32(pR, ro, muR) * sign_(ro);
    }
    if (g2 != 0) {
      ret2 += O32(pR, ro, muR) * ((pN * k) * g_(ro, nu) - k(ro) * pN(nu)) *
              sign_(ro);
    }
    if (g3 != 0) {
      ret3 +=
          O32(pR, ro, muR) * ((k * k) * g_(ro, nu) - k(ro) * k(nu)) * sign_(ro);
    }
  }
  if (spinParity > 0) {
    ret1 = ret1 * ((pR(0) > 0) ? gamma5_ : -gamma5_);
    ret2 = ret2 * gamma5_;
    ret3 = ret3 * gamma5_;
  }
  // cerr << "RNrho: " << g1/(4.*mrho*mrho) << endl;
  // cerr << "g1, mrho: " << g1 << " " << mrho << endl;
  // return ret1;
  return i_ * g1 / (4. * mrho * mrho) * ret1 -
         1. / (8. * mrho * mrho * mrho) * (g2 * ret2 + g3 * ret3);
}

DiracMatrix vertex3hNgamma(double g1, double g2, double g3, halfint spinParity,
                           uint muR, FourVector pR, uint nu, FourVector k) {
  return vertex3hNrho(g1, g2, g3, spinParity, muR, pR, nu, k);
}

DiracMatrix vertexN3hgamma(double g1, double g2, double g3, halfint spinParity,
                           uint muR, FourVector pR, uint nu, FourVector k) {
  return vertexN3hrho(g1, g2, g3, spinParity, muR, pR, nu, k);
}

DiracMatrix vertex5hNpi(double g, halfint spinParity, uint mu, uint nu,
                        FourVector pR, FourVector q) {
  const double mpi = Config::get<double>("pi_pm.mass");
  DiracMatrix ret = g / POW<4>(mpi) *
                    (+g_(mu, nu) * POW<2>(pR * q) + q(mu) * q(nu) * (pR * pR) -
                     (pR(mu) * q(nu) + q(mu) * pR(nu)) * (pR * q)) *
                    gamma_unit;
  if (spinParity > 0) {
    ret = ret * ((pR(0) > 0) ? -gamma5_ : gamma5_);
  }
  return ret;
}

DiracMatrix vertex5hNsi(double g, halfint spinParity, uint mu, uint nu,
                        FourVector pR, FourVector q) {
  const double msi = Config::get<double>("sigma.mass");
  DiracMatrix ret = g / POW<4>(msi) *
                    (+g_(mu, nu) * POW<2>(pR * q) + q(mu) * q(nu) * (pR * pR) -
                     (pR(mu) * q(nu) + q(mu) * pR(nu)) * (pR * q)) *
                    gamma_unit;
  if (spinParity < 0) {
    ret = ret * ((pR(0) > 0) ? -gamma5_ : gamma5_);
  }
  return ret;
}

DiracMatrix vertex5hNrho(double g1, double g2, double g3, halfint spinParity,
                         uint si, uint ta, FourVector pR, uint la,
                         FourVector k) {
  const double mrho = Config::get<double>("rho.mass");
  FourVector pN(-pR - k);
  DiracMatrix ret1(gamma_null);
  double scal2(0);
  double scal3(0);
  for (uint mu : {0, 1, 2, 3}) {
    for (uint nu : {0, 1, 2, 3}) {
      if (g1 != 0) {
        ret1 += O32(k, nu, la) * pN(mu) * O52(pR, mu, nu, si, ta) * sign_(mu) *
                sign_(nu);
      }
      if (g2 != 0) {
        scal2 += (pN(la) * k(nu) - g_(la, nu) * (pN * k)) * pN(mu) *
                 O52(pR, mu, nu, si, ta) * sign_(mu) * sign_(nu);
      }
      if (g3 != 0) {
        scal3 += (k(la) * k(nu) - g_(la, nu) * (k * k)) * pN(mu) *
                 O52(pR, mu, nu, si, ta) * sign_(mu) * sign_(nu);
      }
    }
  }
  DiracMatrix ret2;
  DiracMatrix ret3;
  if (spinParity < 0) {
    ret1 = ret1 * gamma5_;
    ret2 = gamma5_ * scal2;
    ret3 = gamma5_ * scal3;
  } else {
    ret2 = gamma_unit * scal2;
    ret3 = gamma_unit * scal3;
    if (pR(0) > 0) {
      ret2 *= -1.;
      ret3 *= -1.;
    }
  }
  return -i_ / POW<4>(2. * mrho) *
         (+g1 * ret1 + g2 / (2. * mrho) * ret2 + g3 / (2. * mrho) * ret3);
}

DiracMatrix vertex5hNgamma(double g1, double g2, double g3, halfint spinParity,
                           uint si, uint ta, FourVector pR, uint la,
                           FourVector k) {
  return vertex5hNrho(g1, g2, g3, spinParity, si, ta, pR, la, k);
}

/**
   Projectors used in propagators. Terms proportional to the momentum are
   dropped, because the vertices give zero when contracted with the momentum.
*/
DiracMatrix P3h(FourVector p, double m, uint mu, uint nu) {
  assert(isIndex4(mu));
  assert(isIndex4(nu));
  return -g_(mu, nu) * gamma_unit + gamma_(mu) * gamma_(nu) / 3.;
}

DiracMatrix P5h(FourVector p, double m, uint mu, uint nu, uint la, uint ro) {
  assert(isIndex4(mu));
  assert(isIndex4(nu));
  assert(isIndex4(la));
  assert(isIndex4(ro));
  return (-(g_(mu, la) * g_(nu, ro) + g_(mu, ro) * g_(nu, la)) / 2. +
          g_(mu, nu) * g_(la, ro) / 5.) *
             gamma_unit +
         (g_(mu, la) * gamma_(nu) * gamma_(ro) +
          g_(mu, ro) * gamma_(nu) * gamma_(la) +
          g_(nu, la) * gamma_(mu) * gamma_(ro) +
          g_(nu, ro) * gamma_(mu) * gamma_(la)) /
             10.;
}

/**
   Spin 3/2 projector.
*/
DiracMatrix pro3half(FourVector p, double m, uint mu, uint nu) {
  assert(isIndex4(mu));
  assert(isIndex4(nu));
  double m2 = m * m;
  return (-g_(mu, nu) + 2. / (3. * m2) * p(mu) * p(nu)) * gamma_unit +
         gamma_(mu) * gamma_(nu) / 3. +
         (gamma_(mu) * p(nu) - gamma_(nu) * p(mu)) / (3. * m);
}

double _G(FourVector p, double m, uint mu, uint nu) {
  assert(isIndex4(mu));
  assert(isIndex4(nu));
  return -g_(mu, nu) + p(mu) * p(nu) / (m * m);
}

DiracMatrix _T(FourVector p, double m, uint mu, uint nu) {
  const double m2(m * m);
  return -1. / 2. * (gamma_(mu) * gamma_(nu) - gamma_(nu) * gamma_(mu)) +
         1. / (2. * m2) * p(mu) *
             (gamma_(p) * gamma_(nu) - gamma_(nu) * gamma_(p)) -
         1. / (2. * m2) * p(nu) *
             (gamma_(p) * gamma_(mu) - gamma_(mu) * gamma_(p));
}

/**
   Spin 5/2 projector.
*/
DiracMatrix pro5half(FourVector p, double m, uint mu1, uint mu2, uint nu1,
                     uint nu2) {
  assert(isIndex4(mu1));
  assert(isIndex4(mu2));
  assert(isIndex4(nu1));
  assert(isIndex4(nu2));
  return -3. / 10. *
             (_G(p, m, mu1, nu1) * _G(p, m, mu2, nu2) +
              _G(p, m, mu1, nu2) * _G(p, m, mu2, nu1)) *
             gamma_unit +
         1. / 5. * _G(p, m, mu1, mu2) * _G(p, m, nu1, nu2) * gamma_unit +
         1. / 10. *
             (_T(p, m, mu1, nu1) * _G(p, m, mu2, nu2) +
              _T(p, m, mu2, nu2) * _G(p, m, mu1, nu1) +
              _T(p, m, mu1, nu2) * _G(p, m, mu2, nu1) +
              _T(p, m, mu2, nu1) * _G(p, m, mu1, nu2));
}

DiracMatrix fprop(FourVector p, double m) {
  return i_ * (gamma_(p) + m * gamma_unit) / (p * p - m * m);
}

dcomplex sprop(FourVector p, double m) { return i_ / (p * p - m * m); }
