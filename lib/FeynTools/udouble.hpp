#ifndef UDOUBLE_HPP
#define UDOUBLE_HPP

#include <iostream>

class udouble {
public:
  udouble(double value=0, double rel_uncert=0);
  double get_value() const;
  double get_rel_uncert() const;
  double minimum() const;
  double maximum() const;
  friend std::ostream& operator<<(std::ostream&, const udouble&);
private:
  double value;
  double rel_uncert;
};


#endif // UDOUBLE_HPP
