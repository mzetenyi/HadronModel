#include "oldModel.hpp"

#include "Config.hpp"

DiracMatrix vertex3hNpi_old(double g, halfint spinParity, uint muR,
                            FourVector pR, FourVector q) {
  const double mpi = Config::get<double>("pi_pm.mass");
  DiracMatrix ret = -g / mpi * gamma_unit * q(muR);
  return ((spinParity > 0) ? ret : ret * gamma5_);
}

DiracMatrix vertex3hNrho_old(double g, halfint spinParity, uint muR,
                             FourVector pR, uint nu, FourVector k) {
  const double mrho = Config::get<double>("rho.mass");
  DiracMatrix ret =
      -i_ * g / mrho * (gamma_(k) * g_(muR, nu) - k(muR) * gamma_(nu));
  return ((spinParity > 0) ? ret * gamma5_ : ret);
}

DiracMatrix vertex3hNgamma_old(double g, halfint spinParity, uint muR,
                               FourVector pR, uint nu, FourVector k) {
  return vertex3hNrho_old(g, spinParity, muR, pR, nu, k);
}

DiracMatrix P3h_old(FourVector p, double m, uint mu, uint nu) {
  return -g_(mu, nu) * gamma_unit + gamma_(mu) * gamma_(nu) / 3. +
         2. * p(mu) * p(nu) / (3. * m * m) * gamma_unit -
         (p(mu) * gamma_(nu) - p(nu) * gamma_(mu)) / (3. * m);
}