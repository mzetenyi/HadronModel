#ifndef VRANCX_HPP
#define VRANCX_HPP

#include "FeynTools.hpp"
//#include "Vectors.hpp"
using namespace Vectors;
//#include "Spinors.hpp"
using namespace Spinors;
//#include "units.hpp"
using namespace units_GeV;
//#include "Isospin.hpp"

// bool vectorTerm();
// bool tensorTerm();
bool vectorRho();
bool tensorRho();
bool vectorGamma();
bool tensorGamma();
bool vectorOmega();
bool tensorOmega();

double O32(FourVector p, uint mu, uint nu, uint la);
DiracMatrix O32(FourVector p, uint mu, uint la);

double O52(FourVector p, uint mu, uint nu, uint la, uint ro, uint si, uint ta);
double O52(FourVector p, uint mu, uint nu, uint la);

double iso1hh(int Qin, int Qmes);

double isoNNrhopi(int Qin, int Qrho, int Qpi);

dcomplex vertexrhopipi(uint nu, FourVector qp, FourVector qm);

dcomplex vertexgammapipi(uint nu, FourVector qp, FourVector qm);

dcomplex vertexsipipi(FourVector q1, FourVector q2);

DiracMatrix vertexNNpi(FourVector q, int Qin, int Qmes);

DiracMatrix vertexNNsigma();

DiracMatrix vertexNNgamma(uint mu, int Qin, FourVector q);

DiracMatrix vertexNNrho(uint mu, FourVector q, int Qin, int Qmes);

DiracMatrix vertexNNomega(uint mu, FourVector q);

DiracMatrix vertexNNrhopi(uint mu, int Qin, int Qrho, int Qpi);

DiracMatrix vertexNNgammapi(uint mu);

DiracMatrix vertex1hNpi(double g, halfint spinParity, FourVector q);

DiracMatrix vertex1hNsi(double g, halfint spinParity);

DiracMatrix vertex1hNrho(double g1, halfint spinParity, FourVector pR, uint nu, FourVector k);

DiracMatrix vertex1hNgamma(double g1, halfint spinParity, FourVector pR, uint nu, FourVector k);

DiracMatrix vertex3hNpi(double g, halfint spinParity, uint muR, FourVector pR, FourVector q);

DiracMatrix vertex3hNsi(double g, halfint spinParity, uint muR, FourVector pR, FourVector q);

DiracMatrix vertex3hNrho(double g1, double g2, double g3, halfint spinParity, uint muR, FourVector pR, uint nu, FourVector k);
DiracMatrix vertexN3hrho(double g1, double g2, double g3, halfint spinParity, uint muR, FourVector pR, uint nu, FourVector k);

DiracMatrix vertex3hNgamma(double g1, double g2, double g3, halfint spinParity, uint muR, FourVector pR, uint nu, FourVector k);
DiracMatrix vertexN3hgamma(double g1, double g2, double g3, halfint spinParity, uint muR, FourVector pR, uint nu, FourVector k);

DiracMatrix vertex5hNpi(double g, halfint spinParity, uint mu, uint nu, FourVector pR, FourVector q);

DiracMatrix vertex5hNsi(double g, halfint spinParity, uint mu, uint nu, FourVector pR, FourVector q);

DiracMatrix vertex5hNrho(double g1, double g2, double g3, halfint spinParity, uint si, uint ta, FourVector pR, uint la, FourVector k);

DiracMatrix vertex5hNgamma(double g1, double g2, double g3, halfint spinParity, uint si, uint ta, FourVector pR, uint la, FourVector k);

DiracMatrix vertex3h3hpi(double g, 
			 halfint spinParityR, uint muR, FourVector pR, 
			 halfint spinParityD, uint muD, FourVector pD, 
			 FourVector q);

DiracMatrix vertexRNpi(string resonance, FourVector pR, FourVector pN, FourVector q, uint muR1=0, uint muR2=0);

DiracMatrix vertexRNpi(string resonance, FourVector pR, int QR, FourVector pN, int QN, FourVector q, int Qpi, uint muR1=0, uint muR2=0);

DiracMatrix vertexRNgamma(string resonance, FourVector pR, int QR, FourVector pN, FourVector k, uint mu, uint muR1=0, uint muR2=0);

DiracMatrix vertexNRgamma(string resonance, FourVector pR, int QR, FourVector pN, FourVector k, uint mu, uint muR1=0, uint muR2=0);

/** spin-3/2 propagator, without the terms that give 0 contrib. in the case of Vrancx Lagrangians. */
DiracMatrix P3h(FourVector p, double m, uint mu, uint nu);

/** spin-5/2 propagator, without the terms that give 0 contrib. in the case of Vrancx Lagrangians. */
DiracMatrix P5h(FourVector p, double m, uint mu, uint nu, uint la, uint ro);

DiracMatrix pro1half(FourVector p, double m);

DiracMatrix proN(FourVector p);

DiracMatrix pro3half(FourVector p, double m, uint mu, uint nu);

double _G(FourVector p, double m, uint mu, uint nu);

DiracMatrix _T(FourVector p, double m, uint mu, uint nu);

DiracMatrix pro5half(FourVector p, double m, uint mu1, uint mu2, uint nu1, uint nu2);

DiracMatrix fprop(FourVector p, double m);

dcomplex sprop(FourVector p, double m);

#endif // VRANCX_HPP
