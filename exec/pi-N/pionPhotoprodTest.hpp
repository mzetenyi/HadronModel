#ifndef PIONPHOTOPRODTEST
#define PIONPHOTOPRODTEST

#include "Kinema.hpp"
#include "Spinors.hpp"
using namespace Spinors;
#include "Vectors.hpp"
using namespace Vectors;
#include <string>
using namespace std;

double uchCutoff(string resonance, double q);
DiracMatrix P3half(FourVector p, double m, uint mu, uint nu);
DiracMatrix propR(string resonance, FourVector p, uint muR1=0, uint nuR1=0, uint muR2=0, uint nuR2=0);
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
    double threshold() const;
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