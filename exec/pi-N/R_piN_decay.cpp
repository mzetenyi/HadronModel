#include "R_piN_decay.hpp"

#include <string>

#include "Amplitude_R_piN.hpp"
#include "FeynTools.hpp"

using namespace std;
using namespace units_GeV;

double Gamma_R_piN(const std::string& resonance, double srt, int QR) {
  double ret(0);
  for (int QN : {1, 0}) {
    int Qpi = QR - QN;
    if (-1 <= Qpi and Qpi <= 1) {
      Amplitude_R_piN amp(resonance, srt, QR, QN, Qpi);
      DecayWidthSS<Amplitude_R_piN> widthRpiN(amp);
      ret += widthRpiN.GammaTot();
    }
  }
  return ret;
}

void print_Gamma_R_piN(const string& resonance) {
  double mR = getParam<double>(resonance + ".mass");
  double GammaR = getParam<double>(resonance + ".width");
  int QR = getParam<int>("QR", +1);
  vector<double> srtRange =
      getRange("srt", 1 * GeV, mR + 2. * GammaR, 10 * MeV);
  Tabulator tab(9, 14);
  tab.printComment("sqrt(s)", "Gamma");
  tab.printComment("[GeV]", "[GeV]");
  for (double srt : srtRange) {
    tab.printLine(GeV(srt), GeV(Gamma_R_piN(resonance, srt, QR)));
  }
}
