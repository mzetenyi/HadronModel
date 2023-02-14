#include <iomanip>
#include <iostream>

#include "Config.hpp"
#include "units.hpp"
#include "Kinema.hpp"
#include "RandomNumberGenerator.hpp"
#include "Transformations.hpp"
#include "Vectors.hpp"
#include "piN_Npipi.hpp"
#include "piN_Ngammastar.hpp"
using namespace Vectors;
using namespace std;
using namespace units_GeV;

void printMomentum(const FourVector& p) {
  for (int i(0); i < 4; i++) {
    cout << setw(14) << p(i);
  }
  cout << endl;
}

void printGeneratedEvent(double srt) {
  // generate kinematics (M, agles: costhg, phig, costhpi, phipi)
  double mpi = Config::get<double>("pi_pm.mass");
  double mn = Config::get<double>("Nucleon.mass");
  double Mmin = 2. * mpi;
  double Mmax = srt - mn;
  double M = rn(Mmin, Mmax);
  double costhg = rn(-1., 1.);
  double phig = rn(0., 2. * pi_);
  double costhpi = rn(-1., 1.);
  double phipi = rn(0., 2. * pi_);
  // calculate diff. cross section:
  udouble xsec = piN_Npipi_dsigma_dM_dcosthg_dcosthpi_dphipi(srt, M, costhg,
                                                             costhpi, phipi);

  // calculate momenta:
  double pabs = sqrt(lambda(srt * srt, M * M, mn * mn)) / srt;
  double Eg = sqrt(pabs * pabs + M * M);
  double En = sqrt(pabs * pabs + mn * mn);
  double sinthg = sqrt(1. - costhg * costhg);
  double sinthpi = sqrt(1. - costhpi * costhpi);

  ThreeVector e(sinthg * cos(phig), sinthg * sin(phig), costhg);
  FourVector pg_CM(Eg, pabs * e);
  FourVector pn_CM(En, -pabs * e);
  // generate pions in rho rest frame:
  ThreeVector epi(sinthpi * cos(phipi), sinthpi * sin(phipi), costhpi);
  double ppi = momentum(M, mpi, mpi);
  double Epi = sqrt(ppi * ppi + mpi * mpi);
  FourVector ppim_prime(Epi, ppi * epi);
  FourVector ppip_prime(Epi, -ppi * epi);
  // transform to CM frame:
  FourTensor L_CM = TransformationFrom(pg_CM);
  FourVector ppim_CM = L_CM * ppim_prime;
  FourVector ppip_CM = L_CM * ppip_prime;
  // transformation to lab frame:
  double pabs_in = momentum(srt, mn, mpi);
  double En_in = sqrt(mn * mn + pabs_in * pabs_in);
  FourVector pn_in_CM(En_in, 0, 0, -pabs_in);
  FourTensor L_lab = TransformationTo(pn_in_CM);
  FourVector pn = L_lab * pn_CM;
  FourVector pg = L_lab * pg_CM;
  FourVector ppim = L_lab * ppim_CM;
  FourVector ppip = L_lab * ppip_CM;

  // printout:
  cout << "diffsigma: " << mb(xsec.get_value()) << "   +-   " << mb(xsec.get_uncert())
       << endl;
  //cout << setw(10) << "n";
  //printMomentum(pn);
  cout << setw(10) << "pi+";
  printMomentum(ppip);
  cout << setw(10) << "pi-";
  printMomentum(ppim);
}

int main(int argc, char** argv) {
  Config::load(argc, argv);
  double srt(getParam<double>("srt", 1.49));
  // initial printout:
  piN_Ngammastar pNG(srt,0.1);
  vector<string> channels = pNG.getChannels();
  cout << "# events generated for the reaction pi- + p -> n + rho_0 -> n + pi+ + pi-" << endl;
  cout << "# sqrt(s) = " << srt << " GeV" << endl;
  cout << "# M_inv(pi+pi-) is generated uniformly within the kinematical limits" << endl;
  cout << "# virtual rho_0 is generated isotropically in CM frame" << endl;
  cout << "# pi+ and pi- are generated isotropically in rho_0 rest frame" << endl;
  cout << "# differential cross section for each event is given in mb/GeV" << endl;
  cout << "# momentum components are given in GeV in the lab. frame" << endl;
  cout << "# reaction channels included are: ";
  for (string ch : channels) {
    cout << ch << " ";
  }
  cout << endl;
  double phase = getParam<double>("res_rhophase",0);
  cout << "# phase for resonance contributions: " << phase << " deg" << endl;

  int Ngen(getParam<int>("Ngen", 100));
  for (int i(1); i <= Ngen; i++) {
    cout << "Event " << i << endl;
    printGeneratedEvent(srt);
  }
  return 0;
}