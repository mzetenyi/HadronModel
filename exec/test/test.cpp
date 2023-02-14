#include "FeynTools.hpp"

using namespace units_GeV;

int main(int argc, char** argv) {
  double M(1.5 * GeV);
  double m1(900 * MeV);
  double m2(150 * MeV);
  double costh(0.4);
  double phi(1.2);

  Kinema2 kin(M, m1, m2);
  FourVector p(m2, 0, 0, 0);
  FourVector p2 = kin.getp(2, costh, phi);
  FourVector pprime = kin.transformFrom(2, costh, phi, p);
  cout << "p2 = " << p2 << endl;
  cout << "pprime = " << pprime << endl;
  return 0;
}