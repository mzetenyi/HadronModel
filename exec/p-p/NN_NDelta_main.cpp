#include "Amplitude_NN_NDelta.hpp"
#include "CrossSection.hpp"
#include "units.hpp"
using namespace units_GeV;
#include "Tabulator.hpp"

#include <iomanip>
#include <iostream>
using namespace std;

int main(int argc, char** argv) {
  Config::load(argc, argv);
  int QN = getParam<int>("QN", 0);
  int QD = getParam<int>("QD", 2);
  cout << "# NN_NDelta: cross section of p+p -> N + Delta" << endl;
  cout << "\n# Total cross section" << endl;
  Tabulator tab(5,12);
  tab.printComment("sqrt(s)","sigtot");
  tab.printComment("[GeV]","[mb]");
  double srtmin = getParam<double>("srtmin", 2. * GeV);
  double srtmax = getParam<double>("srtmax", 5. * GeV);
  double dsrt = getParam<double>("dsrt", 0.1 * GeV);
  for (double srt(srtmin); srt < srtmax; srt += dsrt) {
    Amplitude_NN_NDelta amplitude(srt, QN, QD);
    CrossSectionSR<Amplitude_NN_NDelta> xsec(amplitude);
    tab.printLine(GeV(srt),mb(xsec.sigmaTot()));
  }
  cout << "\n#parameters:" << endl;
  Config::list(cout, "# ");

  return 0;
}