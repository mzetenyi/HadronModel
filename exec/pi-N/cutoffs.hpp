#ifndef CUTOFFS_HPP
#define CUTOFFS_HPP

#include "halfint.hpp"

#include <string>

/**
   Cutoff factor for R->Npi and R->Nrho vertices. Depends on squared threemom. of the meson in 
   the Resonance rest frame. Equals 1 if q2=0, i.e. NOT on the mass shell!
 */
double cutoffRapp(halfint spin, double q2);

double cutoffZM(double srt, double q2, double mr, double Gamma_Npi, int l);

double cutoff_RNpi(const std::string& resonance, double srt);


#endif // CUTOFFS_HPP
