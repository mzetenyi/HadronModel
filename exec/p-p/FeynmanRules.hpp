#ifndef FEYNMANRULES_HPP
#define FEYNMANRULES_HPP

#include "Config.hpp"
#include "Spinors.hpp"
using namespace Spinors;
#include "Vectors.hpp"
using namespace Vectors;

DiracMatrix vertexNNpi(int Q_Nout, FourVector p_pion, int Q_pion);
double vertexDNpi(int Q_Dout, FourVector p_pion, int Q_pion, int mu);
dcomplex scalarPropagator(FourVector p, double m);

#endif // FEYNMANRULES_HPP