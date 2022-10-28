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
  static const bool tch =
      isSet("tch") or not(isSet("uch") or isSet("interference"));
  static const bool uch =
      isSet("uch") or not(isSet("tch") or isSet("interference"));
  static const bool interference =
      isSet("interference") or not(isSet("tch") or isSet("uch"));
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
  if (isSet("analytic")) {
    double mN = Config::get<double>("mN");
    double mD = Config::get<double>("mD");
    double mN2 = mN * mN;
    double mD2 = mD * mD;
    double mpi2 = mpi * mpi;
    double f_NNpi = Config::get<double>("f_NNpi");
    double f_DNpi = Config::get<double>("f_DNpi");
    double g_NNpi = 2. * mN / mpi * f_NNpi;
    double t = ppi_t * ppi_t;
    double u = ppi_u * ppi_u;
    if (tch) {
      MSQR += POW<4>(formFactor(t)) *
              POW<2>(f_DNpi * g_NNpi / (mpi * (t - mpi2))) * 2. / (3. * mD2) *
              (-t) * (t - POW<2>(mD - mN)) * POW<2>(t - POW<2>(mD + mN));
    }
    if (uch) {
      MSQR += POW<4>(formFactor(u)) *
              POW<2>(f_DNpi * g_NNpi / (mpi * (u - mpi2))) * 2. / (3. * mD2) *
              (-u) * (u - POW<2>(mD - mN)) * POW<2>(u - POW<2>(mD + mN));
    }
    if (interference) {
      MSQR += POW<2>(formFactor(t) * formFactor(u)) *
              POW<2>(f_DNpi * g_NNpi / mpi) * 1. / ((t - mpi2) * (u - mpi2)) *
              1. / (2. * mD2) *
              ((t * u + 2. * (mD2 - mN2) * (t + u) -
                (mD2 - mN2) * (2. * mD2 - mD * mN + 3. * mN2)) *
                   (t * u + mN * (mN + mD) * (mD2 - mN2)) +
               1. / 3. *
                   (t * u + 2. * POW<2>(mD + mN) * (t + u) -
                    POW<2>(mD + mN) * (2. * mD2 + mD * mN + 5. * mN2)) *
                   (t * u - mN * (mD - mN) * (mD2 - mN2)));
    }
  } else {
    for (halfint la1 : {-half, half}) {
      for (halfint la2 : {-half, half}) {
        for (halfint la3 : {-half, half}) {
          for (halfint la4 : {-3 * half, -half, half, 3 * half}) {
            dcomplex helAmp_t(0);
            dcomplex helAmp_u(0);
            for (int mu : {0, 1, 2, 3}) {
              helAmp_t += ubarN(0, la3) * vertexNNpi(QN, -ppi_t, -Qpi_t) *
                          u1(0, la1) * scalarPropagator(ppi_t, mpi) *
                          sign_(mu) * ubarD(mu, la4) *
                          vertexDNpi(QD, ppi_t, Qpi_t, mu) * u2(0, la2);
              helAmp_u += ubarN(0, la3) * vertexNNpi(QN, -ppi_u, -Qpi_u) *
                          u2(0, la2) * scalarPropagator(ppi_u, mpi) *
                          sign_(mu) * ubarD(mu, la4) *
                          vertexDNpi(QD, ppi_u, Qpi_u, mu) * u1(0, la1);
            }
            if (tch) {
              MSQR += real(helAmp_t * conj(helAmp_t));
            }
            if (uch) {
              MSQR += real(helAmp_u * conj(helAmp_u));
            }
            if (interference) {
              MSQR +=
                  real(helAmp_t * conj(helAmp_u) + helAmp_u * conj(helAmp_t));
            }
          }
        }
      }
    }
  }
  return MSQR;
}