#include "FeynTools.hpp"

using namespace units_GeV;

int main(int argc, char** argv) {
  double M(1.5 * GeV);
  double m1(900 * MeV);
  double m2(150 * MeV);
  double costh(.5);
  double phi(1.2);

  Kinema2 kin(M, m1, m2);
  FourVector p(m1, 0, 0, 0);
  FourVector p1 = kin.getp(1, costh, phi);
  FourVector pprime = kin.transformFrom(1, costh, phi, p);
  cout << "p1 = " << p1 << endl;
  cout << "pprime = " << pprime << endl;
  return 0;
}