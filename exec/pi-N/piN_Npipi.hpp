#ifndef PIN_NPIPI_HPP
#define PIN_NPIPI_HPP

#include "udouble.hpp"

udouble piN_Npipi_dsigma_dM(double srt, double M);

udouble piN_Npipi_dsigma_dM_dcosth(double srt, double M, double costh);

udouble piN_Npipi_dsigma_dM_dcosthg_dcosthpi_dphipi(double srt, double M,
                                                    double costhg,
                                                    double costhpi,
                                                    double phipi);

double piN_Npipi_dsigma_dM_from_dilep(double srt, double M);

#endif  // PIN_NPIPI_HPP
