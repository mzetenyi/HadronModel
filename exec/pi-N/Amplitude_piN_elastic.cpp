#include "Amplitude_piN_elastic.hpp"

#include "Config.hpp"
#include "Vectors.hpp"
using namespace Vectors;
#include "Gamma.hpp"
#include "Vrancx.hpp"
#include "wavefunc.hpp"

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
  const std::string res("D1232");
  const double mr = getParam<double>(res + ".mass");
  const double Gr = getParam<double>(res + ".width");
  const int l = getParam<double>(res + ".l");
  const dcomplex BWs = BW(s, mr, Gamma_R(res, s, mr, Gr, l));
  double MSQR(0);
  for (halfint laN_in : {-half, half}) {
    for (halfint laN_out : {-half, half}) {
      dcomplex helAmp(0);
      for (int mu_out : {0, 1, 2, 3}) {
        for (int mu_in : {0, 1, 2, 3}) {
          helAmp +=
              ubarN_out(0, laN_out) *
              vertexRNpi(res, ps, Qs, -pN_out, -QN_out, -ppi_out, -Qpi_out,
                         mu_out) *
              i_ * (gamma_(ps) + gamma_unit * mr) * BWs *
              P3h(ps, mr, mu_out, mu_in) * sign_(mu_out) * sign_(mu_in) *
              vertexRNpi(res, -ps, -Qs, pN_in, QN_in, ppi_in, Qpi_in, mu_in) *
              uN_in(0, laN_in);
        }
        MSQR += real(helAmp * conj(helAmp));
      }
    }
  }
  return MSQR;
}