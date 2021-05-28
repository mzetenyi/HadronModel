#ifndef UNITS_HPP
#define UNITS_HPP

#include <iostream>
#include <map>
#include <string>
#include <cmath>
//#include "debug.hpp"

using namespace std;

/**   Class that represents a unit of measure. Instances correspond to 
   various units. ...
*/
class Unit {
public:
  Unit() : name(""), inMachineUnits(1) {}
  Unit(const string& name, double inMachineUnits);

  double operator*(const Unit& other) const;

  double operator/(const Unit& other) const;

  double operator*(double x) const;

  double operator/(double x) const;

  friend double operator*(double x, const Unit& un);

  friend double operator/(double x, const Unit& un);

  friend double pow(const Unit&, int);

  /**
     Conversion from machine units.
  */
  double operator()(double x) const;

  friend ostream& operator<<(ostream&, const Unit&);
  static double get(const string& name);
  static void printAll(ostream&);
  static double readDimensional(const string&);
private:
  const string name;
  const double inMachineUnits;
  static map<string,double> units;
};

namespace units_GeV {
  extern const double hbarc;
  extern const Unit fm;
  extern const Unit GeV;
  extern const Unit MeV;
  extern const Unit keV;
  extern const Unit TeV;
  extern const Unit fm2;
  extern const Unit fm3;
  extern const Unit meter;
  extern const Unit barn;
  extern const Unit mb;
  extern const Unit mub;
  extern const Unit GeV2;
}

#endif
