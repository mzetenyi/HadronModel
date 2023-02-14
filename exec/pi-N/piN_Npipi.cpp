#include "piN_Npipi.hpp"

#include "Gamma.hpp"
#include "piN_Ngammastar.hpp"
#include "udouble.hpp"

udouble piN_Npipi_dsigma_dM(double srt, double M) {
  piN_Ngammastar DS(srt, M);
  Observable_Xsec_Pionpair_dM OBS(M);
  udouble ret = M / pi_ * DS.expectationValue(OBS);
  return ret;
}

udouble piN_Npipi_dsigma_dM_dcosth(double srt, double M, double costh) {
  piN_Ngammastar DS(srt, M);
  Observable_Xsec_Pionpair_dM_dcosth OBS(M, costh);
  udouble ret = M / pi_ * DS.expectationValue(OBS);
  return ret;
}

udouble piN_Npipi_dsigma_dM_dcosthg_dcosthpi_dphipi(double srt, double M,
                                                    double costhg,
                                                    double costhpi,
                                                    double phipi) {
  piN_Ngammastar DS(srt, M);
  Observable_Xsec_Pionpair_dM_dcosth_dphi OBS(M, costhpi, phipi);
  udouble ret = M / pi_ * DS.expectationValue(costhg, OBS);
  return ret;
}

double piN_Npipi_dsigma_dM_from_dilep(double srt, double M) {
  double grho = Config::get<double>("grho");
  double grho_tilde = Config::get<double>("grho_tilde");
  double mpi = Config::get<double>("pi_pm.mass");
  double mn = Config::get<double>("Nucleon.mass");
  double mrho = Config::get<double>("rho.mass");
  // conv. factor without rhopipi form factor:
  //  double conv_fac =
  //  POW<2>(grho*grho_tilde/(2.*e*e))*POW<3>(sqrt(1.-4.*mpi*mpi/(M*M)));
  // conv factor using mass dependent rho width (which includes form factor):
  double qpi = sqrt(M * M / 4. - mpi * mpi);
  double qpi0 = sqrt(mrho * mrho / 4. - mpi * mpi);
  const double delta(0.2467 * GeV);  // corresponds to  R = 1/delta = 0.8 fm
  double formfac =
      M / mrho * (qpi0 * qpi0 + delta * delta) / (qpi * qpi + delta * delta);
  // double conv_fac =
  // POW<2>(grho*grho_tilde/(2.*e*e))*POW<3>(sqrt(1.-4.*mpi*mpi/(M*M))) *
  // formfac;
  double Gamma_pipi = Gamma_rho(M);
  //  double Gamma_pipi_chk =
  //  POW<2>(grho_tilde)/(48.*pi_)*M*pow(1.-4.*POW<2>(mpi/M),1.5) * formfac;
  //  double Gamma_pipi_chk = POW<2>(grho_tilde)/(6.*pi_)*POW<3>(qpi)/(M*M) *
  //  formfac; cerr << "Gamma_pipi = " << Gamma_pipi << " or " << Gamma_pipi_chk
  //  << endl;
  double Gamma_ee = 1. / (12. * pi_) * POW<2>(e * e / grho) * M;
  double qe = M / 2.;
  //  double conv_fac = qe/qpi * Gamma_pipi/Gamma_ee;
  double conv_fac = Gamma_pipi / Gamma_ee;
  piN_Ngammastar DS(srt, M);
  Observable_Xsec_Dilepton_dM OBS(M);
  double ret = DS.expectationValue(OBS).get_value() * conv_fac;
  return ret;
}
