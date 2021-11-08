#ifndef OLDMODEL_HPP
#define OLDMODEL_HPP

#include "Vectors.hpp"
using namespace Vectors;
#include "Spinors.hpp"
using namespace Spinors;

DiracMatrix vertex3hNpi_old(double g, halfint spinParity, uint muR, FourVector pR, FourVector q);

DiracMatrix vertex3hNrho_old(double g1, double g2, double g3, halfint spinParity, uint muR, FourVector pR, uint nu, FourVector k);

DiracMatrix vertex3hNgamma_old(double g1, double g2, double g3, halfint spinParity, uint muR, FourVector pR, uint nu, FourVector k);


#endif // OLDMODEL_HPP
