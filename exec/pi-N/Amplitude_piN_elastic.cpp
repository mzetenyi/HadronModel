#include "Amplitude_piN_elastic.hpp"

#include "Config.hpp"
#include "Vectors.hpp"
using namespace Vectors;
#include "FormFactors.hpp"
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
  if (isSet("analytic")) return MSQR_analytic(costh);
  return MSQR_numeric(costh);
}

double Amplitude_piN_elastic::MSQR_numeric(double costh) {
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
  int npol = 2;
  return ret / npol;
}

double Amplitude_piN_elastic::MSQR_analytic(double costh) {
  FourVector pN_in = kinIn.p1();
  FourVector ppi_in = kinIn.p2();
  FourVector pN_out = kinOut.p1(costh);
  FourVector ppi_out = kinOut.p2(costh);
  int Qs = QN_in + Qpi_in;
  FourVector ps = pN_in + ppi_in;
  ubar_ ubarN_out(half, pN_out);
  u_ uN_in(half, pN_in);
  const double s = srt * srt;
  string res = "N1440";

  const double mpi = getParam<double>("pi_pm.mass");
  const double mR = getParam<double>(res + ".mass");
  const double Gr = getParam<double>(res + ".width");
  const int l = getParam<double>(res + ".l");
  const dcomplex BWs = BW(s, mR, Gamma_R(res, s, mR, Gr, l));
  const double BW2 = real(BWs * conj(BWs));
  double g = Config::get<double>(res + ".g0");
  double pNi_pNo = pN_in * pN_out;
  double pNi_ppii = pN_in * ppi_in;
  double pNo_ppii = pN_out * ppi_in;

  double ret =

      +pNi_pNo *
          (-4. * POW<6>(mpi) + 4. * POW<2>(mR) * POW<4>(mpi) -
           8. * mN * mR * POW<4>(mpi) - 12. * POW<2>(mN) * POW<4>(mpi) -
           16. * POW<2>(mN) * POW<2>(mR) * POW<2>(mpi) -
           32. * POW<3>(mN) * mR * POW<2>(mpi) -
           16. * POW<4>(mN) * POW<2>(mpi) - 16. * pNi_ppii * POW<4>(mpi) +
           16. * pNi_ppii * POW<2>(mR) * POW<2>(mpi) -
           16. * pNi_ppii * POW<2>(mN) * POW<2>(mpi) -
           16. * POW<2>(pNi_ppii) * POW<2>(mpi) +
           16. * POW<2>(pNi_ppii) * POW<2>(mR) +
           32. * POW<2>(pNi_ppii) * mN * mR +
           16. * POW<2>(pNi_ppii) * POW<2>(mN))

      + 4. * POW<2>(mN) * POW<6>(mpi) +
      4. * POW<2>(mN) * POW<2>(mR) * POW<4>(mpi) +
      24. * POW<3>(mN) * mR * POW<4>(mpi) + 20. * POW<4>(mN) * POW<4>(mpi) +
      16. * POW<4>(mN) * POW<2>(mR) * POW<2>(mpi) +
      32. * POW<5>(mN) * mR * POW<2>(mpi) + 16. * POW<6>(mN) * POW<2>(mpi) +
      8. * pNi_ppii * pNo_ppii * POW<4>(mpi) -
      8. * pNi_ppii * mN * mR * POW<4>(mpi) +
      8. * pNi_ppii * POW<2>(mN) * POW<4>(mpi) +
      32. * pNi_ppii * POW<3>(mN) * mR * POW<2>(mpi) +
      32. * pNi_ppii * POW<4>(mN) * POW<2>(mpi) +
      32. * POW<2>(pNi_ppii) * pNo_ppii * POW<2>(mpi) +
      32. * POW<2>(pNi_ppii) * pNo_ppii * mN * mR +
      32. * POW<2>(pNi_ppii) * pNo_ppii * POW<2>(mN) -
      32. * POW<2>(pNi_ppii) * mN * mR * POW<2>(mpi) -
      16. * POW<2>(pNi_ppii) * POW<2>(mN) * POW<2>(mpi) -
      16. * POW<2>(pNi_ppii) * POW<2>(mN) * POW<2>(mR) -
      32. * POW<2>(pNi_ppii) * POW<3>(mN) * mR -
      16. * POW<2>(pNi_ppii) * POW<4>(mN) + 32. * POW<3>(pNi_ppii) * pNo_ppii -
      32. * POW<3>(pNi_ppii) * mN * mR - 32. * POW<3>(pNi_ppii) * POW<2>(mN) -
      8. * pNo_ppii * mN * mR * POW<4>(mpi) -
      8. * pNo_ppii * POW<2>(mN) * POW<4>(mpi) -
      16. * pNo_ppii * POW<2>(mN) * POW<2>(mR) * POW<2>(mpi) -
      32. * pNo_ppii * POW<3>(mN) * mR * POW<2>(mpi) -
      16. * pNo_ppii * POW<4>(mN) * POW<2>(mpi);

  ret *= POW<4>(g / mpi);
  double FF = formfactorRNpi(res, ps * ps);
  double isofac_in = (Qpi_in == 0 ? 1. : sqrt(2.));
  double isofac_out = (Qpi_out == 0 ? 1. : sqrt(2.));
  int npol(2);
  return ret * BW2 * POW<2>(FF * FF * isofac_in * isofac_out) / npol;
}
