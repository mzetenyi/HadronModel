#ifndef PIONPHOTOPRODTEST
#define PIONPHOTOPRODTEST

#include "Kinema.hpp"
#include "Spinors.hpp"
using namespace Spinors;

DiracMatrix vertexRNpi(FourVector pR, FourVector pN, FourVector q);
DiracMatrix vertexRNgamma(FourVector pR, FourVector pN, FourVector k, uint mu);
DiracMatrix pro1half(FourVector p, double m);
double N1440width(double srt);
dcomplex BreitWigner(FourVector p, double m, double Gamma);
DiracMatrix propR(FourVector p);
DiracMatrix proN(FourVector p);

class pionPhotoprodTest {
public:
    pionPhotoprodTest(double srt);
    double MSQRraw_analytic(double costh);
    double MSQRraw_numeric(double costh);
    double diffsig_analytic(double costh);
    double diffsig_numeric(double costh);
    double sigtot_analytic();
    double sigtot_numeric();
private:
    double srt;
    double mR;
    double mN;
    double mpi;
    Kinema2 KINin;
    Kinema2 KINout;
};

#endif // PIONPHOTOPRODTEST