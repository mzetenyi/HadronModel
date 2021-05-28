#ifndef GAMMA_HPP
#define GAMMA_HPP

#include <string>

#include "utils.hpp"

double Gamma_rho(double m); ///< Mass dependence of total width = width of rho -> pi + pi

double Gamma_R(const std::string& res, double s, double m0, double G0, int l);

dcomplex BW(double s, double M, double Gamma);

class GammaResonance {
public:
  GammaResonance(const std::string& resonance);
  double operator()(double srt) const;
private:
  std::string res;
  double mass;
  double Gamma0;
  int l;
};

#endif // GAMMA_HPP
