#ifndef OLDMODEL_HPP
#define OLDMODEL_HPP

#include "Vectors.hpp"
using namespace Vectors;
#include "Spinors.hpp"
using namespace Spinors;

DiracMatrix vertex3hNpi_old(double g, halfint spinParity, uint muR, FourVector pR, FourVector q);

DiracMatrix vertex3hNrho_old(double g, halfint spinParity, uint muR, FourVector pR, uint nu, FourVector k);

DiracMatrix vertex3hNgamma_old(double g, halfint spinParity, uint muR, FourVector pR, uint nu, FourVector k);

/** spin-3/2 propagator, INCLUDING the terms that give 0 contrib. in the case of Vrancx Lagrangians. */
DiracMatrix P3h_old(FourVector p, double m, uint mu, uint nu);

#endif // OLDMODEL_HPP
