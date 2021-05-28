#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include <iostream>
#include <string>
#include <vector>

class Histogram1d {
public:
  Histogram1d(double min, double max, uint n);
  void add(double val, double weight=1.);
  void print(std::ostream&, std::string label_val, std::string label_weight, uint width_val=12, uint width_weight=12);
private:
  double min;
  double max;
  uint n;
  std::vector<double> content;
};


#endif // HISTOGRAM_HPP
