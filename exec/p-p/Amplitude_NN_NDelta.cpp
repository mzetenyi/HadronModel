#include "Amplitude_NN_NDelta.hpp"

#include "Config.hpp"
#include "Vectors.hpp"
using namespace Vectors;
#include "FeynmanRules.hpp"
#include "wavefunc.hpp"

Amplitude_NN_NDelta::Amplitude_NN_NDelta(double srt, int QN, int QD)
    : srt(srt), QN(QN), QD(QD) {}

Kinema2 Amplitude_NN_NDelta::getInputKinematics() {
  double mN = Config::get<double>("mN");
  return Kinema2(srt, mN, mN);
}

Kinema2 Amplitude_NN_NDelta::getOutputKinematics(double m) {
  double mN = Config::get<double>("mN");
  return Kinema2(srt, mN, m);
}

double Amplitude_NN_NDelta::width(double m) {
  double Gamma = Config::get<double>("GammaD");
  return Gamma;
}

double Amplitude_NN_NDelta::spectralFunction(double m) {
  double mD = Config::get<double>("mD");
  double Gamma = width(m);
  return 1. / pi_ * mD * Gamma /
         (POW<2>(m * m - mD * mD) + mD * mD * Gamma * Gamma);
}

double Amplitude_NN_NDelta::MSQR(double m, double costh) {
  double mpi = Config::get<double>("mpi");
  Kinema2 kinIn = getInputKinematics();
  Kinema2 kinOut = getOutputKinematics(m);
  FourVector p1 = kinIn.p1();
  FourVector p2 = kinIn.p2();
  FourVector p3 = kinOut.p1(costh);
  FourVector p4 = kinOut.p2(costh);
  FourVector ppi_t = p1 - p3;
  FourVector ppi_u = p2 - p3;
  int Q1 = 1;
  int Q2 = 1;
  int Q3 = QN;
  int Q4 = QD;
  int Qpi_t = Q1 - Q3;
  int Qpi_u = Q2 - Q3;
  ubar_ ubarN(half, p3);
  ubar_ ubarD(3 * half, p4);
  u_ u1(half, p1);
  u_ u2(half, p2);
  double MSQR(0);
  for (halfint la1 : {-half, half}) {
    for (halfint la2 : {-half, half}) {
      for (halfint la3 : {-half, half}) {
        for (halfint la4 : {-3 * half, -half, half, 3 * half}) {
          dcomplex helAmp(0);
          for (int mu : {0, 1, 2, 3}) {
            helAmp += ubarN(0, la3) * vertexNNpi(QN, -ppi_t, -Qpi_t) *
                      u1(0, la1) * scalarPropagator(ppi_t, mpi) * sign_(mu) *
                      ubarD(mu, la4) * vertexDNpi(QD, ppi_t, Qpi_t, mu) *
                      u2(0, la2);
            helAmp += ubarN(0, la3) * vertexNNpi(QN, -ppi_u, -Qpi_u) *
                      u2(0, la2) * scalarPropagator(ppi_u, mpi) * sign_(mu) *
                      ubarD(mu, la4) * vertexDNpi(QD, ppi_u, Qpi_u, mu) *
                      u1(0, la1);
          }
          MSQR += real(helAmp * conj(helAmp));
        }
      }
    }
  }
  return MSQR;
}