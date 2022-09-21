#include "piN_elastic.hpp"

#include "Amplitude_piN_elastic.hpp"
#include "FeynTools.hpp"

using namespace units_GeV;

void print_piN_elastic_sigtot() {
  vector<double> srtRange = getRange("srt", 1 * GeV, 2 * GeV, 10 * MeV);
  int QN_in = getParam<int>("QN_in", +1);
  int Qpi_in = getParam<int>("Qpi_in", -1);
  int QN_out = getParam<int>("QN_out", QN_in);
  int Qpi_out = getParam<int>("Qpi_out", QN_in + Qpi_in - QN_out);
  Tabulator tab(9, 14);
  tab.printComment("sqrt(s)", "sigma_tot");
  tab.printComment("[GeV]", "[mb]");
  for (double srt : srtRange) {
    Amplitude_piN_elastic amp(srt, QN_in, Qpi_in, QN_out, Qpi_out);
    CrossSectionSS<Amplitude_piN_elastic> xsec(amp);
    double sigtot = xsec.sigmaTot();
    tab.printLine(GeV(srt), mb(sigtot));
  }
}

void print_piN_elastic_dsig() {
  vector<double> costhRange = getRange("costh", -1, 1.0001, 0.1);
  int QN_in = getParam<int>("QN_in", +1);
  int Qpi_in = getParam<int>("Qpi_in", -1);
  int QN_out = getParam<int>("QN_out", QN_in);
  int Qpi_out = getParam<int>("Qpi_out", QN_in + Qpi_in - QN_out);
  double srt = getParam<double>("srt", 1.5 * GeV);
  Tabulator tab(13, 14);
  tab.printComment("cos(theta)", "dsigma/dcosth");
  tab.printComment("", "[mb]");
  Amplitude_piN_elastic amp(srt, QN_in, Qpi_in, QN_out, Qpi_out);
  CrossSectionSS<Amplitude_piN_elastic> xsec(amp);
  for (double costh : costhRange) {
    double dsig = xsec.dsigma_dcosth(costh);
    tab.printLine(costh, mb(dsig));
  }
}