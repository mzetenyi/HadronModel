#include "FormFactors.hpp"
#include "Config.hpp"
#include "Kinema.hpp"

using namespace std;

double formfactorRNpi(string resonance, double m2) {
  if (Config::exists("noRNpiFF")) return 1;
  double m = (m2 > 0) ? sqrt(m2) : 0;
  double mN = Config::get<double>("Nucleon.mass");
  double mpi = Config::get<double>("pi_pm.mass");
  double mR = Config::get<double>(resonance + ".mass");
  double Gamma = Config::get<double>(resonance + ".width");
  int l = Config::get<double>(resonance + ".l");
  double delta2 = pow(mR - mN - mpi, 2) + Gamma * Gamma / 4.;
  double q0 = momentum(mR, mN, mpi);
  double q = (m > mN + mpi) ? momentum(m, mN, mpi) : 0;
  // return sqrt(mR / m) *
  //       pow((q0 * q0 + delta2) / (q * q + delta2), (l + 1.) / 2.);
  // double ret= sqrt(mR / m) *
  //       pow((q0 * q0 + delta2) / (q * q + delta2), (l + 1.) / 2.);
  double ret = pow((q0 * q0 + delta2) / (q * q + delta2),
                   (l + 1.) / 2.);  // correspond to Deniz
  if (isinf(ret)) {
    cerr << "inf in formfactorRNpi" << endl;
    PR(sqrt(mR / m));
    PR(pow((q0 * q0 + delta2) / (q * q + delta2), (l + 1.) / 2.));
    exit(0);
  }
  // PR(ret);
  return ret;
}
