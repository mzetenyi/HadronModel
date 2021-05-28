#include "units.hpp"

#include <iostream>
#include <map>
#include <string>
#include <cmath>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
using boost::algorithm::split;
using boost::algorithm::is_any_of;
using boost::algorithm::token_compress_on;
using boost::algorithm::trim;

using namespace std;

Unit::Unit(const string& name, double inMachineUnits) :
  name(name),
  inMachineUnits(inMachineUnits) {
  units[name] = inMachineUnits;
}

double Unit::operator*(const Unit& other) const {
  return inMachineUnits * other.inMachineUnits;
}

double Unit::operator/(const Unit& other) const {
  return inMachineUnits / other.inMachineUnits;
}

double Unit::operator*(double x) const {
  return inMachineUnits*x;
}

double Unit::operator/(double x) const {
  return inMachineUnits/x;
}

double operator*(double x, const Unit& un) {
  return x*un.inMachineUnits;
}

double operator/(double x, const Unit& un) {
  return x/un.inMachineUnits;
}

double pow(const Unit& un, int i) {
  return pow(un.inMachineUnits,i);
}

/**
   Conversion from machine units.
*/
double Unit::operator()(double x) const {
  return x/inMachineUnits;
}


ostream& operator<<(ostream& out, const Unit& un) {
  out << "Unit: " << un.name << " = " << un.inMachineUnits << endl;
  return out;
}

double Unit::get(const string& name) {
  return units.at(name);
}

void Unit::printAll(ostream& out) {
  for (auto x : units) {                                                                                   
    out << x.first << " = " << x.second << endl;
  }                                                                                                       
  out << endl;
}

double Unit::readDimensional(const string& str) {
  vector<string> tokens;
  split(tokens,str,is_any_of("*"),token_compress_on);
  double val = stod(tokens[0]);
  if (tokens.size()>1) {
    trim(tokens[1]);
    double un = get(tokens[1]);
    val *= un;
  }
  return val;
}



map<string,double> Unit::units;

namespace units_GeV {
  const double hbarc = 0.197327; // GeV*fm
  const Unit fm("fm",1./hbarc);
  const Unit GeV("GeV",1.);
  const Unit MeV("MeV",GeV/1000.);
  const Unit keV("keV",MeV/1000.);
  const Unit TeV("TeV",1000.*GeV);
  const Unit fm2("fm2",fm*fm);
  const Unit fm3("fm3",fm2*fm);
  const Unit meter("meter",1e15*fm);
  const Unit barn("barn",1e-28*meter*meter);
  const Unit mb("mb",1e-3*barn);
  const Unit mub("mub",1e-6*barn);
  const Unit GeV2("GeV2",GeV*GeV);
}

