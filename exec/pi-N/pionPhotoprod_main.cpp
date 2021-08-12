#include <iomanip>
#include <iostream>

#include "Config.hpp"
#include "pionPhotoprod.hpp"
#include "udouble.hpp"
#include "units.hpp"
using namespace units_GeV;

using namespace std;

void print_xsec(ostream& out) {
  static const double mn = Config::get<double>("Nucleon.mass");
  static const double mpi = Config::get<double>("pi_pm.mass");
  double srtmin(1.08 * GeV);
  double srtmax(1.8 * GeV);
  double dsrt(0.01 * GeV);

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
  out << "#" << setw(12) << "srt [GeV]"
      << " " << setw(15) << "sig_tot [mub]"
      << " " << setw(15) << "rel_error"
      << " " << setw(15) << "min"
      << " " << setw(15) << "max" << endl;
  for (double srt(srtmin); srt < srtmax; srt += dsrt) {
    udouble sig = pionPhotoprod_sigtot(srt);
    out << setw(13) << srt << " " << setw(15) << mub(sig.get_value()) << " "
        << setw(15) << sig.get_rel_uncert() << " " << setw(15)
        << max(mub(sig.minimum()), 0.) << " " << setw(15) << mub(sig.maximum())
        << endl;
  }
}

int main(int argc, char** argv) {
  Config::load(argc, argv);

  print_xsec(cout);

  if (not Config::exists("-n")) {
    cout << "#===============================================================\n"
         << endl;
    cout << "# Parameters:" << endl;
    Config::list(cout, "# ");
    cout << "#===============================================================\n"
         << endl;
    cout << "#" << endl;
    ;
  }

  return 1;
}
