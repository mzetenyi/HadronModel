#include "FeynmanRules.hpp"
#include "Isospin.hpp"

DiracMatrix vertexNNpi(int Q_Nout, FourVector p_pion, int Q_pion) {
  double f_NNpi = Config::get<double>("f_NNpi");
  double m_pi = Config::get<double>("m_pi");
  return f_NNpi / m_pi * gamma5_ * gamma_(p_pion) *
         isospin_1h1h1(Q_Nout, Q_Nout - Q_pion, Q_pion);
}

double vertexDNpi(int Q_Dout, FourVector p_pion, int Q_pion, int mu) {
  double f_DNpi = Config::get<double>("f_DNpi");
  double m_pi = Config::get<double>("m_pi");
  return f_DNpi / m_pi * p_pion(mu) *
         isospin_3h1h1(Q_Dout, Q_Dout - Q_pion, Q_pion);
}

dcomplex scalarPropagator(FourVector p, double m) {
  return i_ / (p * p - m * m);
}