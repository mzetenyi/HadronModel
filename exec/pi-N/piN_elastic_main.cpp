#include "FeynTools.hpp"
#include "piN_elastic.hpp"

int main(int argc, char** argv) {
  Config::load(argc,argv);
  if (isSet("sigtot")) {
    print_piN_elastic_sigtot();
  }
  if (isSet("dsig")) {
    cerr << "calculating dsig" << endl;
    print_piN_elastic_dsig();
  }
  return 0;
}