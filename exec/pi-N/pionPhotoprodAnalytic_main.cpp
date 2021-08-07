#include "pionPhotoprodAnalytic.hpp"
#include <iostream>
#include <iomanip>
#include "udouble.hpp"
#include "Config.hpp"
#include "units.hpp"
using namespace units_GeV;

using namespace std;

void print_xsec(ostream& out) {
  static const double mn = Config::get<double>("Nucleon.mass");
  static const double mpi = Config::get<double>("pi_pm.mass");
  double srtmin(1.08*GeV);
  double srtmax(1.8*GeV);
  double dsrt(0.01*GeV);

  if (Config::exists("srtmin")) {
    srtmin = Config::get<double>("srtmin");
  }
  if (Config::exists("srtmax")) {
    srtmax = Config::get<double>("srtmax");
  }
  if (Config::exists("dsrt")) {
    dsrt = Config::get<double>("dsrt");
  }
  out << "# cross section of pion photoproduction" << endl;
  out << "# using the analytic calculation for s-channel N(1440) contribution" << endl;
  out << "#" << setw(12) << "srt [GeV]"
      << " " << setw(15) << "sig_tot [mub]"
      << endl;
  for (double srt(srtmin); srt<srtmax; srt+=dsrt) {
    pionPhotoprodAnalytic PPA(srt);
    double sig = PPA.sigtot();
    out << setw(13) << srt
        << " " << setw(15) << mub(sig)
        << endl;
  }
}


int main(int argc, char** argv) {

  Config::load(argc,argv);
  
  print_xsec(cout);

  if (not Config::exists("-n")) {
    cout << "# ======================================================================\n" << endl;
    cout << "# Parameters:" << endl;
    Config::list(cout,"# ");
    cout << "# ======================================================================\n" << endl;
    cout << "#" << endl;;
  }
  
  return 1;
}
