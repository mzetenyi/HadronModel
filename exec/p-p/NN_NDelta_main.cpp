#include "Amplitude_NN_NDelta.hpp"
#include "CrossSection.hpp"
#include "units.hpp"
using namespace units_GeV;
#include <iomanip>
#include <iostream>

#include "Tabulator.hpp"
using namespace std;

int main(int argc, char** argv) {
  Config::load(argc, argv);
  int QN = getParam<int>("QN", 0);
  int QD = getParam<int>("QD", 2);
  double srt = getParam<double>("srt", 2.5);
  cout << "# NN_NDelta: cross section of p+p -> N + Delta" << endl;
  cout << "\n# Diff cross section dsigma/dmD" << endl;
  cout << "# sqrt(s) = " << srt << " GeV" << endl;
  Tabulator tab(5, 12);
  tab.printComment("mD", "dsigma/dmD");
  tab.printComment("[GeV]", "[mb]");
  double mpi = getParam<double>("mpi");
  double mN = getParam<double>("mN");
  double mDmin = mN + mpi;
  double mDmax = srt - mN;
  double mDstart = 1. * GeV;
  double dmD = getParam<double>("dmD", 0.1 * GeV);
  for (double mD(mDstart); mD < mDmax; mD += dmD) {
    if (mD > mDmin) {
      Amplitude_NN_NDelta amplitude(srt, QN, QD);
      CrossSectionSR<Amplitude_NN_NDelta> xsec(amplitude);
      tab.printLine(GeV(mD), mb(xsec.dsigma_dm(mD)));
    }
  }
  cout << "\n#parameters:" << endl;
  Config::list(cout, "# ");

  return 0;
}