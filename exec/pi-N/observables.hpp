#ifndef OBSERVABLES_HPP
#define OBSERVABLES_HPP

#include "utils.hpp"
#include "MultiArray.hpp"

class Observable_Xsec_Dilepton_dM_dcosth_dphi {
public:
  Observable_Xsec_Dilepton_dM_dcosth_dphi() = delete;
  Observable_Xsec_Dilepton_dM_dcosth_dphi(double M, double costh, double phi);
  dcomplex operator()(halfint lap, halfint la) const;
private:
  double M;
  double costh;
  double phi;
  MultiArray<dcomplex> obs; 
};

class Observable_Xsec_Dilepton_dM_dcosth {
public:
  Observable_Xsec_Dilepton_dM_dcosth() = delete;
  Observable_Xsec_Dilepton_dM_dcosth(double M, double costh);
  dcomplex operator()(halfint lap, halfint la) const;
private:
  double M;
  double costh;
  MultiArray<dcomplex> obs; 
};

class Observable_Xsec_Dilepton_dM {
public:
  Observable_Xsec_Dilepton_dM() = delete;
  Observable_Xsec_Dilepton_dM(double M);
  dcomplex operator()(halfint lap, halfint la) const;
private:
  double M;
  MultiArray<dcomplex> obs; 
};

class Observable_Xsec_Pionpair_dcosth_dphi {
public:
  Observable_Xsec_Pionpair_dcosth_dphi() = delete;
  Observable_Xsec_Pionpair_dcosth_dphi(double M, double costh, double phi);
  dcomplex operator()(halfint lap, halfint la) const;
private:
  double M;
  double costh;
  double phi;
  MultiArray<dcomplex> obs; 
};

class Observable_Xsec_Pionpair_dcosth {
public:
  Observable_Xsec_Pionpair_dcosth() = delete;
  Observable_Xsec_Pionpair_dcosth(double M, double costh);
  dcomplex operator()(halfint lap, halfint la) const;
private:
  double M;
  double costh;
  MultiArray<dcomplex> obs; 
};

class Observable_Xsec_Pionpair {
public:
  Observable_Xsec_Pionpair() = delete;
  Observable_Xsec_Pionpair(double M);
  dcomplex operator()(halfint lap, halfint la) const;
private:
  double M;
  MultiArray<dcomplex> obs; 
};


#endif // OBSERVABLES_HPP
