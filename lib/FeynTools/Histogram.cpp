#include "Histogram.hpp"
#include <iostream>
#include <iomanip>
#include <cassert>
//#include "utils.hpp"

using namespace std;

Histogram1d::Histogram1d(double min, double max, uint n) : 
min(min),
  max(max),
  n(n),
  content(n,0.) {
}

void Histogram1d::add(double val, double weight) {
  assert(min<val);
  assert(val<max);
  // cerr << "Histogram1d::add ---------------------------------------" << endl;
  // cerr << "val = " << val << endl;
  // cerr << "weight = " << weight << endl;
  // PR(min); PR(max); PR((val-min)/(max-min)); PR(int((val-min)/(max-min)));
  content[int(n*(val-min)/(max-min))] += weight;
}

void Histogram1d::print(ostream& out, string label_val, string label_weight,
			uint width_val, uint width_weight) {
  out << setw(2*width_val+1) << label_val << " " << setw(width_weight) << label_weight << endl;
  double dval = (max-min)/n;
  for (uint i(0); i<n; i++) {
    out << setw(width_val) << min+i*dval << " " 
	<< setw(width_val) << min+(i+1)*dval << " " 
	<< setw(width_weight) << content[i] << endl;
  }
}
