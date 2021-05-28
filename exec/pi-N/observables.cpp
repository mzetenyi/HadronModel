#include "observables.hpp"

#include "Config.hpp"
#include "units.hpp"

using namespace units_GeV;

Observable_Xsec_Dilepton_dcosth_dphi::Observable_Xsec_Dilepton_dcosth_dphi(double M, double costh, double phi) :
  M(M),
  costh(costh),
  phi(phi),
  obs(idx_s1,idx_s1) {
  double sinth = sqrt(1.-costh*costh);
  double fac = M*M;
  obs( _1, _1) =   fac * (1. + costh*costh);
  obs( _1, _0) = - fac * sqrt(2.)*costh*sinth*exp(-i_*phi);
  obs( _1,-_1) =   fac * sinth*sinth*exp(-2.*i_*phi);
  obs( _0, _1) = - fac * sqrt(2.)*costh*sinth*exp(i_*phi);
  obs( _0, _0) =   fac * 2.*sinth*sinth;
  obs( _0,-_1) =   fac * sqrt(2.)*costh*sinth*exp(-i_*phi);
  obs(-_1, _1) =   fac * sinth*sinth*exp(2.*i_*phi);
  obs(-_1, _0) =   fac * sqrt(2.)*costh*sinth*exp(i_*phi);
  obs(-_1,-_1) =   fac * (1. + costh*costh);
}

dcomplex Observable_Xsec_Dilepton_dcosth_dphi::operator()(halfint lap, halfint la) const {
  return obs(lap,la);
}

Observable_Xsec_Dilepton_dcosth::Observable_Xsec_Dilepton_dcosth(double M, double costh) :
  M(M),
  costh(costh),
  obs(idx_s1,idx_s1) {
  double sinth = sqrt(1.-costh*costh);
  double fac = 2.*pi_ * M*M;
  obs( _1, _1) = fac * (1. + costh*costh);
  obs( _1, _0) = 0;
  obs( _1,-_1) = 0;
  obs( _0, _1) = 0;
  obs( _0, _0) = fac * 2.*sinth*sinth;
  obs( _0,-_1) = 0;
  obs(-_1, _1) = 0;
  obs(-_1, _0) = 0;
  obs(-_1,-_1) = fac * (1. + costh*costh);
}

dcomplex Observable_Xsec_Dilepton_dcosth::operator()(halfint lap, halfint la) const {
  return obs(lap,la);
}

Observable_Xsec_Dilepton::Observable_Xsec_Dilepton(double M) :
  M(M),
  obs(idx_s1,idx_s1) {
  double fac = 2.*pi_ * 8./3. * M*M;
  obs( _1, _1) = fac;
  obs( _1, _0) = 0;
  obs( _1,-_1) = 0;
  obs( _0, _1) = 0;
  obs( _0, _0) = fac;
  obs( _0,-_1) = 0;
  obs(-_1, _1) = 0;
  obs(-_1, _0) = 0;
  obs(-_1,-_1) = fac;
}

dcomplex Observable_Xsec_Dilepton::operator()(halfint lap, halfint la) const {
  return obs(lap,la);
}

Observable_Xsec_Pionpair_dcosth_dphi::Observable_Xsec_Pionpair_dcosth_dphi(double M, double costh, double phi) :
  M(M),
  costh(costh),
  phi(phi),
  obs(idx_s1,idx_s1) {
  double sinth = sqrt(1.-costh*costh);
  double mpi = Config::get<double>("pi_pm.mass");
  double grho = Config::get<double>("grho");
  double grho_tilde = Config::get<double>("grho_tilde");
  double corrfac = POW<2>(grho/(e*e)); // remove the rho-gamma transition vertex from piN_Ngammastar - must use the rho (and not rho2) contrib!
  corrfac *= sqrt(1. - 4.*mpi*mpi/(M*M)); // correct for difference in phase-space factor abs(q_pi) or abs (q_e)
  const double delta(0.2467*GeV);  // corresponds to  R = 1/delta = 0.8 fm
  double mrho = Config::get<double>("rho.mass");
  double qpi = sqrt(M*M/4. - mpi*mpi);
  double qpi0 = sqrt(mrho*mrho/4. - mpi*mpi);
  double formfac = M/mrho * (qpi0*qpi0 + delta*delta)/(qpi*qpi + delta*delta);
  double fac = (M < 2.*mpi) ? 0 : 2.*POW<2>(grho_tilde)*(M*M/4. - mpi*mpi) * formfac * corrfac; // 2*pabs * formfac * corrfac
  obs( _1, _1) =   fac * (1. - costh*costh);
  obs( _1, _0) =   fac * sqrt(2.)*costh*sinth*exp(-i_*phi);
  obs( _1,-_1) = - fac * sinth*sinth*exp(-2.*i_*phi);
  obs( _0, _1) =   fac * sqrt(2.)*costh*sinth*exp(i_*phi);
  obs( _0, _0) =   fac * 2.*costh*costh;
  obs( _0,-_1) = - fac * sqrt(2.)*costh*sinth*exp(-i_*phi);
  obs(-_1, _1) = - fac * sinth*sinth*exp(2.*i_*phi);
  obs(-_1, _0) = - fac * sqrt(2.)*costh*sinth*exp(i_*phi);
  obs(-_1,-_1) =   fac * (1. - costh*costh);
}

dcomplex Observable_Xsec_Pionpair_dcosth_dphi::operator()(halfint lap, halfint la) const {
  return obs(lap,la);
}

Observable_Xsec_Pionpair_dcosth::Observable_Xsec_Pionpair_dcosth(double M, double costh) :
  M(M),
  costh(costh),
  obs(idx_s1,idx_s1) {
  double sinth = sqrt(1.-costh*costh);
  double mpi = Config::get<double>("pi_pm.mass");
  double grho = Config::get<double>("grho");
  double grho_tilde = Config::get<double>("grho_tilde");
  double corrfac = POW<2>(grho/(e*e)); // remove the rho-gamma transition vertex from piN_Ngammastar - must use the rho (and not rho2) contrib!
  corrfac *= sqrt(1. - 4.*mpi*mpi/(M*M)); // correct for difference in phase-space factor abs(q_pi) or abs (q_e)
  const double delta(0.2467*GeV);  // corresponds to  R = 1/delta = 0.8 fm
  double mrho = Config::get<double>("rho.mass");
  double qpi = sqrt(M*M/4. - mpi*mpi);
  double qpi0 = sqrt(mrho*mrho/4. - mpi*mpi);
  double formfac = M/mrho * (qpi0*qpi0 + delta*delta)/(qpi*qpi + delta*delta);
  double fac = (M < 2.*mpi) ? 0 : 2.*pi_ * 2.*POW<2>(grho_tilde)*(M*M/4. - mpi*mpi) * formfac * corrfac;
  obs( _1, _1) = fac * (1. - costh*costh);
  obs( _1, _0) = 0;
  obs( _1,-_1) = 0;
  obs( _0, _1) = 0;
  obs( _0, _0) = fac * 2.*costh*costh;
  obs( _0,-_1) = 0;
  obs(-_1, _1) = 0;
  obs(-_1, _0) = 0;
  obs(-_1,-_1) = fac * (1. - costh*costh);
}

dcomplex Observable_Xsec_Pionpair_dcosth::operator()(halfint lap, halfint la) const {
  return obs(lap,la);
}

Observable_Xsec_Pionpair::Observable_Xsec_Pionpair(double M) :
  M(M),
  obs(idx_s1,idx_s1) {
  double mpi = Config::get<double>("pi_pm.mass");
  double grho = Config::get<double>("grho");
  double grho_tilde = Config::get<double>("grho_tilde");
  double corrfac = POW<2>(grho/(e*e)); // remove the rho-gamma transition vertex from piN_Ngammastar - must use the rho (and not rho2) contrib!
  corrfac *= sqrt(1. - 4.*mpi*mpi/(M*M)); // correct for difference in phase-space factor abs(q_pi) or abs (q_e)
  const double delta(0.2467*GeV);  // corresponds to  R = 1/delta = 0.8 fm
  double mrho = Config::get<double>("rho.mass");
  double qpi = sqrt(M*M/4. - mpi*mpi);
  double qpi0 = sqrt(mrho*mrho/4. - mpi*mpi);
  double formfac = M/mrho * (qpi0*qpi0 + delta*delta)/(qpi*qpi + delta*delta);
  double fac = (M < 2.*mpi) ? 0 : 2.*pi_ * 4./3. * 2.*POW<2>(grho_tilde)*(M*M/4. - mpi*mpi) * formfac * corrfac;
  obs( _1, _1) = fac;
  obs( _1, _0) = 0;
  obs( _1,-_1) = 0;
  obs( _0, _1) = 0;
  obs( _0, _0) = fac;
  obs( _0,-_1) = 0;
  obs(-_1, _1) = 0;
  obs(-_1, _0) = 0;
  obs(-_1,-_1) = fac;
}

dcomplex Observable_Xsec_Pionpair::operator()(halfint lap, halfint la) const {
  return obs(lap,la);
}

