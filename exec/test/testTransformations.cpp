#include "FeynTools.hpp"

int main(int argc, char** argv) {
  double m(1.2);
  FourVector p_rest(m, 0, 0, 0);
  ThreeVector P(0.1, 0.2, 0.3);
  double pabs = sqrt(P * P);
  double E = sqrt(m * m + P * P);
  FourVector p_boost(E, P);
  FourTensor L = TransformationTo(p_boost);
  FourTensor Linv = TransformationFrom(p_boost);
  FourVector p_trans = L * p_boost;
  FourVector p_invtrans = Linv * p_rest;
  cout << "p_rest = " << p_rest << endl;
  cout << "p_boost = " << p_boost << endl;
  cout << "p_trans = " << p_trans << endl;
  cout << "p_invtrans = " << p_invtrans << endl;
  cout << "L * Linv = " << L * Linv << endl;
  return 0;
}