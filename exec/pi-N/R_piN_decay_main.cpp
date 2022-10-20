#include <string>

#include "FeynTools.hpp"
#include "R_piN_decay.hpp"

using namespace std;

int main(int argc, char** argv) {
  Config::load(argc, argv);
  string resonance = getParam<string>("resonance");
  print_Gamma_R_piN(resonance);
  return 0;
}