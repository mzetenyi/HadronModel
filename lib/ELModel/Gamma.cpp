#include <cmath>
#include <iostream>
#include <string>
#include "Gamma.hpp"
#include "Config.hpp"
#include "units.hpp"
using namespace units_GeV;


double Gamma_rho(double m) { ///< Mass dependence of total width = width of rho -> pi + pi
  const double delta(0.2467*GeV);  // corresponds to  R = 1/delta = 0.8 fm 
  static const double mpi_pm = Config::get<double>("pi_pm.mass");
  static const double mpi = mpi_pm;
  static const double m0 = Config::get<double>("rho.mass");
  static const double Gamma0 = Config::get<double>("rho.width");
  double x = m*m - 4.*mpi*mpi;
  if (x<=0) { return 0; }
  double q = sqrt(x)/2.;  // pion momentum
  double q0 = sqrt(m0*m0 - 4.*mpi*mpi)/2.; // pion momentum at the pole
  int rhospec = 1;
  if (Config::exists("rhospec")) { rhospec = Config::get<int>("rhospec"); }
  if (rhospec==0) { return Gamma0 * m0/m * pow(q/q0,3) * (q0*q0 + delta*delta)/(q*q + delta*delta); } // Manley, also used by BoGa, with R = 1/delta = 0.8 fm
  if (rhospec==1) { return Gamma0 * pow(q/q0,3); } // version used with Enrico
  if (rhospec==2) { return Gamma0 * m0/m * pow(q/q0,3); } // Nucl.Phys. A632 (1998) 109-127 (Giessen)
  if (rhospec==3) { return Gamma0 * pow(q*m/(q0*m0),3); } // version used in Deniz' code
  std::cerr << "Unknown width parametrization of rho width, rhospec=" << rhospec << std::endl;
  exit(0);
  return 1; // never get here
}

double Gamma_R(const std::string& res, double s, double m0, double G0, int l) {
  if (Config::exists("fixedRwidth")) return G0;
  static const double mpi_pm = Config::get<double>("pi_pm.mass");
  //  static const double mpi_0 = Config::get<double>("pi_0.mass");
  //static const double mpi = (2.*mpi_pm + mpi_0)/3.;
  static const double mpi = mpi_pm;
  static const double mn = Config::get<double>("Nucleon.mass");
  if ((s>0) and (sqrt(s)>mn+mpi)) {
    double q  = sqrt(lambda(s,    mn*mn,mpi*mpi))/(2.*sqrt(s));
    double q0 = sqrt(lambda(m0*m0,mn*mn,mpi*mpi))/(2.*m0);
    double d2;
    if (res == "D1232") {
      d2 = POW<2>(0.3);
    } else if (res == "N1535") {
      d2 = POW<2>(0.5);
    } else {
      d2 = POW<2>(m0-mn-mpi) + G0*G0/4.;
    }
    return G0 * m0/sqrt(s) * pow(q/q0,2*l+1) * pow((q0*q0+d2)/(q*q+d2),l+1);
  }
  return 0;
}

dcomplex BW(double s, double M, double Gamma) {
  if (s>0) {
    return 1./(s - M*M + i_*sqrt(s)*Gamma);
    //return 1./(s - M*M + i_*M*Gamma);
  }
  return 1./(s - M*M);
}

GammaResonance::GammaResonance(const std::string& resonance) :
  res(resonance),
  mass(Config::get<double>(resonance+".mass")),
  Gamma0(Config::get<double>(resonance+".width")),
  l(Config::get<int>(resonance+".l")) {
}

double GammaResonance::operator()(double srt) const {
  return Gamma_R(res,srt*srt, mass, Gamma0, l);
}
