#ifndef PIONPHOTOPRODTEST
#define PIONPHOTOPRODTEST

#include "Kinema.hpp"
#include "Spinors.hpp"
using namespace Spinors;
#include "Vectors.hpp"
using namespace Vectors;
#include <string>
using namespace std;

double formfactorRNpi(string resonance, double m);
DiracMatrix vertexRNpi(FourVector pR, FourVector pN, FourVector q);
DiracMatrix vertexRNgamma(FourVector pR, FourVector pN, FourVector k, uint mu);
DiracMatrix vertexRNpi(string resonance, FourVector pR, FourVector pN, FourVector q, uint muR1=0, uint muR2=0);
DiracMatrix vertexRNgamma(string resonance, FourVector pR, FourVector pN, FourVector k, uint mu, uint muR1=0, uint muR2=0);
DiracMatrix pro1half(FourVector p, double m);
double N1440width(double srt);
double resonanceWidth(string resonance, double m);
dcomplex BreitWigner(FourVector p, double m, double Gamma);
DiracMatrix propR(FourVector p);
DiracMatrix propR(string resonance, FourVector p, uint muR1=0, uint nuR1=0, uint muR2=0, uint nuR2=0);
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