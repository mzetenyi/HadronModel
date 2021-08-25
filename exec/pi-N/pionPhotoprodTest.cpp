#include "pionPhotoprodTest.hpp"

#include "Config.hpp"
#include "units.hpp"
using namespace units_GeV;

#include "Spinors.hpp"
#include "utils.hpp"
using namespace Spinors;

#include "MultiArray.hpp"

DiracMatrix vertexRNpi(FourVector pR, FourVector pN, FourVector q) {
  return gamma5_ * gamma_(q);
}

DiracMatrix vertexRNgamma(FourVector pR, FourVector pN, FourVector k, uint mu) {
  return gamma_(k) * gamma_(mu) - gamma_(mu) * gamma_(k);
}

DiracMatrix propR(FourVector p) {
  double mR(Config::get<double>("N1440.mass"));
  return gamma_(p) + mR * gamma_unit;
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
  PR(mR);
  PR(mN);
  PR(mpi);
}

double pionPhotoprodTest::MSQRraw_numeric(double costh) {
  FourVector pi = KINin.p1(1);
  FourVector k = KINin.p2(1);
  FourVector pf = KINout.p1(costh);
  FourVector q = KINout.p2(costh);
  FourVector p = pi + k;

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

  MultiArray<DiracMatrix> T(idx_lor);
  for (uint mu(0); mu < 4; mu++) {
    T(mu) = vertexRNpi(p, pf, q) * propR(p) * vertexRNgamma(p, pi, k, mu);
  }
  dcomplex MSQR(0);
  for (uint mu(0); mu < 4; mu++) {
    MSQR += trace(proN(pf) * T(mu) * proN(pi) * adj(T(mu)) * sign_(mu));
  }
  return real(MSQR);
}

double pionPhotoprodTest::MSQRraw_analytic(double costh) {
  FourVector pi = KINin.p1(1);
  FourVector k = KINin.p2(1);
  FourVector pf = KINout.p1(costh);
  FourVector q = KINout.p2(costh);
  FourVector p = pi + k;

  PR(pi);
  PR(k);
  PR(pf);
  PR(q);

  double pi_pf = pi * pf;
  double pf_q = pf * q;
  double pi_q = pi * q;
  double pi_k = pi * k;
  double k_p = k * p;
  double p_pf = p * pf;
  double p_q = p * q;
  double k_pf = k * pf;
  double k_q = k * q;

  double MSQR = +128 * pi_k * k_p * p_pf * POW<2>(mpi) -
                256 * pi_k * k_p * p_q * pf_q +
                128 * pi_k * k_p * mR * mN * POW<2>(mpi) -
                64 * pi_k * k_pf * POW<2>(mpi) * POW<2>(srt) +
                64 * pi_k * k_pf * POW<2>(mR) * POW<2>(mpi) +
                128 * pi_k * k_q * pf_q * POW<2>(srt) -
                128 * pi_k * k_q * pf_q * POW<2>(mR);

  return MSQR;
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
  return 1;
}