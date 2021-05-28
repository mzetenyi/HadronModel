#include "piN_Ndilep.hpp"
#include <iostream>
#include <iomanip>
#include "udouble.hpp"
#include "Config.hpp"
#include "units.hpp"
using namespace units_GeV;

using namespace std;

void print_invarian_mass_spectrum(ostream& out, double srt) {
  static const double mn = Config::get<double>("Nucleon.mass");
  double Mmin(0.);
  double Mmax(srt-mn);
  double dM(0.01*GeV);

  if (Config::exists("Mmin")) {
    Mmin = Config::get<double>("Mmin");
  }
  if (Config::exists("Mmax")) {
    Mmax = Config::get<double>("Mmax");
  }
  if (Config::exists("dM")) {
    dM = Config::get<double>("dM");
  }
  out << "# dilepton invariant mass spectrum" << endl;
  out << "#" << setw(12) << "M [GeV]"
      << " " << setw(15) << "dsig_dM [mub/GeV]"
      << " " << setw(15) << "rel_error"
      << " " << setw(15) << "min"
      << " " << setw(15) << "max"
      << endl;
  for (double M(Mmin); M<Mmax; M+=dM) {
    udouble dsig = piN_Ndilep_dsigma_dM(srt,M);
    out << setw(13) << M
        << " " << setw(15) << mub(dsig.get_value())
        << " " << setw(15) << dsig.get_rel_uncert()
        << " " << setw(15) << max(mub(dsig.minimum()),0.)
        << " " << setw(15) << mub(dsig.maximum())
        << endl;
  }
}


int main(int argc, char** argv) {

  Config::load(argc,argv);
  double srt = Config::get<double>("srt");
  
  print_invarian_mass_spectrum(cout,srt);

  if (not Config::exists("-n")) {
    cout << "# ======================================================================\n" << endl;
    cout << "# Parameters:" << endl;
    Config::list(cout,"# ");
    cout << "# ======================================================================\n" << endl;
    cout << "#" << endl;;
  }
  
  return 1;
}
