#include "Resonance.hpp"

#include "FeynTools.hpp"

using namespace std;

vector<string> getResonances() {
  vector<string> allres = {"D1232", "N1440", "N1520",
                           "N1535", "N1675", "N1680"};
  vector<string> ret;
  for (string res : allres) {
    if (isSet(res)) ret.push_back(res);
  }
  return ret;
}

vector<int> getIndexRange(string resonance, int i) {
  halfint spin = getParam<halfint>(resonance + ".spin");
  if (i <= itrunc(spin)) return {0, 1, 2, 3};
  return {0};
}