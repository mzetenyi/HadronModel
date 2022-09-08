#ifndef RESDEC_HPP
#define RESDEC_HPP

#include <string>
#include <functional>

double Gsigma_pipi(double g, double mm);

double GammaR_Npi(std::string res, double mr);

double dGR_Nrho_dm(std::string res, double mr, double mm); // mm: mass of rho

double GammaR_Nrho(std::string res, double mr);

double GammaR_Ngamma(std::string res, double mr, std::string VMDversion="VMD1");

double Gamma_sigma(double m); ///< Mass dependence of total width = width of sigma -> pi + pi

double dGR_Nsigma_dm(std::string res, double mr, double mm); // mm: mass of sigma

double GammaR_Nsigma(std::string res, double mr);

double dGR_Dpi_dmdcosth(std::string res, double mr, double mm, double costh); // mr: res mass, mm: running Delta mass

double dGR_Dpi_dm(std::string res, double mr, double mm);

double GammaR_Dpi(std::string res, double mr);

double GammaR_tot(std::string res, double mr);

class bindResonance {
public:
  bindResonance(std::function<double(std::string,double)> func, std::string resonance);
  double operator()(double) const;
private:
  std::function<double(std::string,double)> func;
  std::string resonance;
  double scaleFactor;
  static double getScaleFactor(std::function<double(std::string,double)> func, std::string resonance);
};

#endif // RESDEC_HPP
