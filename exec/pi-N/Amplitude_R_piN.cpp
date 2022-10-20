#include "Amplitude_R_piN.hpp"

#include "Config.hpp"
#include "Resonance.hpp"
#include "Vrancx.hpp"
#include "utils.hpp"
#include "wavefunc.hpp"

using namespace std;

Amplitude_R_piN::Amplitude_R_piN(string resonance, double srt, int QR, int QN,
                                 int Qpi)
    : resonance(resonance),
      M(getParam<double>(resonance + ".mass")),
      srt(srt),
      QR(QR),
      QN(QN),
      Qpi(Qpi),
      mN(getParam<double>("Nucleon.mass")),
      mpi(getParam<double>(Qpi == 0 ? "pi_0.mass" : "pi_pm.mass")),
      kin(srt, mN, mpi) {}

Kinema2 Amplitude_R_piN::getKinematics() {
  return kin;
}

double Amplitude_R_piN::MSQR() {
  vector<halfint> laRrange = getHelicityRange(resonance);
  vector<int> mu1range = getIndexRange(resonance, 1);
  vector<int> mu2range = getIndexRange(resonance, 2);
  FourVector pR = kin.P();
  FourVector pN = kin.p1();
  FourVector ppi = kin.p2();
  ubar_ ubar_N(half, pN);
  halfint spinR = getParam<halfint>(resonance + ".spin");
  u_ u_R(spinR, pR);
  double ret(0);
  for (halfint laR : laRrange) {
    for (halfint laN : {-half, half}) {
      dcomplex helamp(0);
      for (int mu1 : mu1range) {
        for (int mu2 : mu2range) {
          helamp +=
              ubar_N(laN) *
              vertexRNpi(resonance, pR, QR, -pN, -QN, -ppi, -Qpi, mu1, mu2) *
              u_R(laR, mu1, mu2) * sign_(mu1) * sign_(mu2);
        }
      }
      ret += real(helamp * conj(helamp));
    }
  }
  return ret;
}
