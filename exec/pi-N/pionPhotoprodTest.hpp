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
double uchCutoff(string resonance, double q);
DiracMatrix vertexRNpi(string resonance, FourVector pR, FourVector pN, FourVector q, uint muR1=0, uint muR2=0);
DiracMatrix vertexRNpi(string resonance, FourVector pR, int QR, FourVector pN, int QN, FourVector q, int Qpi, uint muR1=0, uint muR2=0);
DiracMatrix vertexRNgamma(string resonance, FourVector pR, int QR, FourVector pN, FourVector k, uint mu, uint muR1=0, uint muR2=0);
DiracMatrix pro1half(FourVector p, double m);
double resonanceWidth(string resonance, double m);
dcomplex BreitWigner(string resonance, double srt);
DiracMatrix propR(string resonance, FourVector p, uint muR1=0, uint nuR1=0, uint muR2=0, uint nuR2=0);
DiracMatrix proN(FourVector p);
double widthRNpi(string resonance, double M);
double widthRNpi(string resonance);
double widthRNpi(string resonance, int QR, int QN, int Qpi, double M);
double widthRNpi(string resonance, int QR, double M);
double widthRNpi(string resonance, int QR, int QN, int Qpi);
double widthRNpi(string resonance, int QR);
double widthRNgamma(string resonance, int QR, double M);
double widthRNgamma(string resonance, int QR);

class pionPhotoprodTest {
public:
    pionPhotoprodTest(double srt, int Qpi, int Qf);
    //double MSQRraw_analytic(double costh);
    double MSQR_numeric(double costh);
    //double diffsig_analytic(double costh);
    double diffsig_numeric(double costh);
    //double sigtot_analytic();
    double sigtot_numeric();
private:
    double srt;
    int Qpi;
    int Qf;
    int Qi;
    double mR;
    double mN;
    double mp;
    double mn;
    double mNi;
    double mNf;
    double mpipm;
    double mpi0;
    double mpi;
    Kinema2 KINin;
    Kinema2 KINout;
};

#endif // PIONPHOTOPRODTEST