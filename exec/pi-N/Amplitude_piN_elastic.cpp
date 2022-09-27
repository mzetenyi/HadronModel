#include "Amplitude_piN_elastic.hpp"

#include "Config.hpp"
#include "Vectors.hpp"
using namespace Vectors;
#include "Gamma.hpp"
#include "Resonance.hpp"
#include "Vrancx.hpp"
#include "wavefunc.hpp"

using namespace std;

Amplitude_piN_elastic::Amplitude_piN_elastic(double srt, int QN_in, int Qpi_in,
                                             int QN_out, int Qpi_out)
    : srt(srt),
      QN_in(QN_in),
      Qpi_in(Qpi_in),
      QN_out(QN_out),
      Qpi_out(Qpi_out),
      mN(getParam<double>("Nucleon.mass")),
      mpi_in(getParam<double>(Qpi_in == 0 ? "pi_0.mass" : "pi_pm.mass")),
      mpi_out(getParam<double>(Qpi_out == 0 ? "pi_0.mass" : "pi_pm.mass")),
      kinIn(Kinema2(srt, mN, mpi_in)),
      kinOut(Kinema2(srt, mN, mpi_out)) {}

Kinema2 Amplitude_piN_elastic::getInputKinematics() { return kinIn; }

Kinema2 Amplitude_piN_elastic::getOutputKinematics() { return kinOut; }

double Amplitude_piN_elastic::MSQR(double costh) {
  FourVector pN_in = kinIn.p1();
  FourVector ppi_in = kinIn.p2();
  FourVector pN_out = kinOut.p1(costh);
  FourVector ppi_out = kinOut.p2(costh);
  int Qs = QN_in + Qpi_in;
  FourVector ps = pN_in + ppi_in;
  ubar_ ubarN_out(half, pN_out);
  u_ uN_in(half, pN_in);
  const double s = srt * srt;
  vector<string> resonances = getResonances();
  double ret(0);
  for (halfint laN_in : {-half, half}) {
    for (halfint laN_out : {-half, half}) {
      dcomplex helAmp(0);
      for (string res : resonances) {
        const double mr = getParam<double>(res + ".mass");
        const double Gr = getParam<double>(res + ".width");
        const int l = getParam<double>(res + ".l");
        const dcomplex BWs = BW(s, mr, Gamma_R(res, s, mr, Gr, l));
        vector<int> range1 = getIndexRange(res, 1);
        vector<int> range2 = getIndexRange(res, 2);
        for (int mu1_out : range1) {
          for (int mu2_out : range2) {
            for (int mu1_in : range1) {
              for (int mu2_in : range2) {
                helAmp += ubarN_out(0, laN_out) *
                          vertexRNpi(res, ps, Qs, -pN_out, -QN_out, -ppi_out,
                                     -Qpi_out, mu1_out, mu2_out) *
                          i_ *
                          resonanceProjector(res, ps, mu1_out, mu2_out, mu1_in,
                                             mu2_in) *
                          BWs * sign_(mu1_out) * sign_(mu2_out) *
                          sign_(mu1_in) * sign_(mu2_in) *
                          vertexRNpi(res, -ps, -Qs, pN_in, QN_in, ppi_in,
                                     Qpi_in, mu1_in, mu2_in) *
                          uN_in(0, laN_in);
              }
            }
          }
        }
      }
      ret += real(helAmp * conj(helAmp));
    }
  }
  return ret;
}