#ifndef BORNTERMS_HPP
#define BORNTERMS_HPP

#include "utils.hpp"
#include "MultiArray.hpp"
#include "Vectors.hpp"
using namespace Vectors;
#include "Spinors.hpp"
using namespace Spinors;


class BornTerms {
public:
  BornTerms(FourVector p1, int Q1, FourVector q, int Qpi, FourVector k);
  double Fs;
  double Fu;
  double Ft;
  double Fhat;
  MultiArray<DiracMatrix> M_;
  MultiArray<DiracMatrix> Tgamma_s;
  MultiArray<DiracMatrix> Tgamma_u;
  MultiArray<DiracMatrix> Tgamma_t;
  MultiArray<DiracMatrix> Tgamma_c;
  MultiArray<DiracMatrix> Tgamma_corr;
  MultiArray<DiracMatrix> Tgamma_all;
  MultiArray<DiracMatrix> Trho_all;
  dcomplex Agamma[6];
  dcomplex Arho[6];
};



#endif // BORNTERMS_HPP
