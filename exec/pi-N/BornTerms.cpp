#include "BornTerms.hpp"
#include "Gamma.hpp"
#include "Config.hpp"


// All momenta and charges are incoming!
BornTerms::BornTerms(FourVector p1, int Q1, FourVector q, int Qpi, FourVector k) :
  Fs(1.),
  Fu(1.),
  Ft(1.),
  Fhat(1.),
  M_(idx(_1,_5),idx_lor),
  Tgamma_s(idx_lor),
  Tgamma_u(idx_lor),
  Tgamma_t(idx_lor),
  Tgamma_c(idx_lor),
  Tgamma_corr(idx_lor),
  Tgamma_all(idx_lor),
  Trho_all(idx_lor) {
  FourVector p2 = - p1 - q - k;
  int Q2 = - Q1 - Qpi;
  double ISO = -1. + (Q1*Q1 + Q2*Q2) + sqrt(2.)*Qpi*Qpi;
  static const double f_NNpi = Config::get<double>("f_NNpi");
  static const double mpi = Config::get<double>("pi_pm.mass");
  static const double mn = Config::get<double>("Nucleon.mass");
  static const double mrho = Config::get<double>("rho.mass");
  static const double kan = Config::get<double>("ka_nngamma");
  static const double kap = Config::get<double>("ka_ppgamma");
  static const double kar = Config::get<double>("ka_NNrho");
  // static const double kan = 0;
  // static const double kap = 0;
  // static const double kar = 0;
  static const double ka1 = (Q1==0) ? kan : kap;
  static const double ka2 = (Q2==0) ? kan : kap;
  static const double alpha(1/137.036);
  static const double e(sqrt(4.*pi_*alpha));
  static const double fac = e*f_NNpi/mpi * ISO;
  if (Config::exists("cutoffBorn")) {
    //    cerr << "CUTOFF" << endl;
    double LBorn = 0.63; // from Zetenyi, Wolf PRC 2012
    if (Config::exists("LBorn")) {
      LBorn = Config::get<double>("LBorn");
    }

    Fs = 1./(1 + POW<2>(2.*(p2*k)+k*k)/POW<4>(LBorn));
    Fu = 1./(1 + POW<2>(2.*(p1*k)+k*k)/POW<4>(LBorn));
    Ft = 1./(1 + POW<2>(2.*(q*k)+k*k)/POW<4>(LBorn));
    Fhat = Fs + Fu + Ft - Fs*Fu - Fs*Ft - Fu*Ft + Fs*Fu*Ft;
  }
  for (halfint mu(0); mu<4; mu++) {
    M_(_1,mu) = gamma5_*(gamma_(mu)*gamma_(k) - k(mu)*gamma_unit);
    M_(_2,mu) = gamma5_/2.*((2.*p1(mu)+k(mu))*(2.*(q*k)+k*k) - (2.*q(mu)+k(mu))*(2.*(p1*k)+k*k));
    M_(_3,mu) = gamma5_/2.*((2.*p2(mu)+k(mu))*(2.*(q*k)+k*k) - (2.*q(mu)+k(mu))*(2.*(p2*k)+k*k));
    M_(_4,mu) = gamma5_/2.*(gamma_(mu)*(2.*(p1*k)+k*k) - (2.*p1(mu)+k(mu))*gamma_(k));
    M_(_5,mu) = gamma5_/2.*(gamma_(mu)*(2.*(p2*k)+k*k) - (2.*p2(mu)+k(mu))*gamma_(k));
  }
  Agamma[1] = fac * (2.*mn*((Q2+ka2)*Fs/(2.*(p2*k)+k*k) - (Q1-ka1)*Fu/(2.*(p1*k)+k*k)) +
                     (ka1*Fu + ka2*Fs)/(2.*mn));
  Agamma[2] = fac * 4.*mn * Q1*Fhat/(2.*(p1*k)+k*k)/(2.*(q*k)+k*k);
  Agamma[3] = fac * 4.*mn * Q2*Fhat/(2.*(p2*k)+k*k)/(2.*(q*k)+k*k);
  Agamma[4] = -fac * 2.*ka1*Fu/(2.*(p1*k)+k*k);
  Agamma[5] = fac * 2.*ka2*Fs/(2.*(p2*k)+k*k);
  for (halfint mu(0); mu<4; mu++) {
    Tgamma_s(mu) = Fs * fac * (Q2*gamma5_*gamma_(mu)
                          + Q2/(2*(p2*k)+k*k) *
                          2.*mn*gamma5_*(gamma_(mu)*gamma_(k) + 2.*p2(mu)*gamma_unit));
    Tgamma_u(mu) = Fu * fac * (Q1*gamma5_*gamma_(mu)
                          + Q1/(2*(p1*k)+k*k) *
                          2.*mn*gamma5_*(gamma_(k)*gamma_(mu) + 2.*p1(mu)*gamma_unit));
    Tgamma_t(mu) = Ft * fac * Qpi/(2*(q*k)+k*k) * 2.*mn*gamma5_*(2.*q(mu)+k(mu));
    Tgamma_c(mu) = Ft * fac * Qpi*gamma5_*gamma_(mu);
    Tgamma_corr(mu) = - fac * 2.*mn*gamma5_ * (Q1*(Fu-Fhat)*(2.*p1(mu)+k(mu))/(2.*(p1*k)+k*k)
                                             + Q2*(Fs-Fhat)*(2.*p2(mu)+k(mu))/(2.*(p2*k)+k*k)
                                             + Qpi*(Ft-Fhat)*(2.*q(mu)+k(mu))/(2.*(q*k)+k*k))
      - fac * gamma5_*gamma_(mu) * (Q2*(Fs-Fhat) + Q1*(Fu-Fhat) + Qpi*(Ft-Fhat));
    Tgamma_all(mu) = gamma_null;
    for (int j(1); j<=5; j++) {
      Tgamma_all(mu) += Agamma[j] * M_(_1*j,mu);
    }
  }
  const double grho_tilde = Config::get<double>("grho_tilde");
  const double grho = Config::get<double>("grho");
  const double M2 = k*k;
  const double M = sqrt(M2);
  const dcomplex F_rho = grho_tilde/grho * (M2)/(M2 - mrho*mrho + i_*sqrt(M2)*Gamma_rho(M));
  Arho[1] = fac*F_rho * (2.*mn*((Q2+0.5)*(1.-kar)*Fs/(2.*(p2*k)+k*k) - (Q1-0.5)*(1.+kar)*Fu/(2.*(p1*k)+k*k)) -
                           kar/(2.*mn) * (Fs*(Q2+0.5) - Fu*(Q1-0.5)));
  Arho[2] = fac*F_rho * 4.*mn * (Q1-0.5)*Fhat/(2.*(p1*k)+k*k)/(2.*(q*k)+k*k);
  Arho[3] = fac*F_rho * 4.*mn * (Q2+0.5)*Fhat/(2.*(p2*k)+k*k)/(2.*(q*k)+k*k);
  Arho[4] = -fac*F_rho * 2.*kar*Fu*(Q1-0.5)/(2.*(p1*k)+k*k);
  Arho[5] = -fac*F_rho * 2.*kar*Fs*(Q2+0.5)/(2.*(p2*k)+k*k);
  for (halfint mu(0); mu<4; mu++) {
    Trho_all(mu) = gamma_null;
    for (int j(1); j<=5; j++) {
      Trho_all(mu) += Arho[j] * M_(_1*j,mu);
    }
  }  
}
