#ifndef PIN_NDILEP_HPP
#define PIN_NDILEP_HPP

#include "udouble.hpp"

udouble piN_Ndilep_dsigma_dM(double srt, double M);
udouble piN_Ndilep_dsigma_dM_dcosth(double srt, double M, double costh);

/**
   Integrated cross section above Mthr.
 */
double piN_Ndilep_sigma_highM(double srt, double Mthr);

#endif // PIN_NDILEP_HPP
