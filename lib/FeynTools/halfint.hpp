#ifndef HALFINT_HPP
#define HALFINT_HPP

#include <iostream>

using std::ostream;
using std::istream;

#include "utils.hpp"

/**
   A representation of half-integers, like 1/2, 2/2, 3/2, etc.
*/
class halfint {
public:
  /** Default constructor, initialize to 0. */
  halfint();
  /** Convert an int to halfint. The value of the resulting halfint will be x. */
  explicit halfint(int x);
  /** copy constructor */
  halfint(const halfint&);
  halfint& operator=(const halfint&);
  halfint& operator=(int);
  /** Conversion to double. */
  //operator double() const { return _2x/2.; }
  const halfint& operator+() const;
  const halfint operator-() const;
  const halfint& operator++();
  const halfint operator++(int);
  const halfint& operator--();
  const halfint operator--(int);
  const halfint& operator+=(const halfint&);
  const halfint& operator-=(const halfint&);
  const halfint& operator*=(int);
  const double operator*(double);
  const double operator/(double);
  const dcomplex operator*(dcomplex);
  const dcomplex operator/(dcomplex);
  friend const halfint operator+(const halfint&, const halfint&);
  friend const halfint operator+(const halfint&, int);
  friend const halfint operator+(int, const halfint&);
  friend double operator+(const halfint&, double);
  friend double operator+(double, const halfint&);
  friend dcomplex operator+(const halfint&, dcomplex);
  friend dcomplex operator+(dcomplex, const halfint&);
  friend const halfint operator-(const halfint&, const halfint&);
  friend const halfint operator-(const halfint&, int);
  friend const halfint operator-(int, const halfint&);
  friend double operator-(const halfint&, double);
  friend double operator-(double, const halfint&);
  friend dcomplex operator-(const halfint&, dcomplex);
  friend dcomplex operator-(dcomplex, const halfint&);
  friend const double operator*(const halfint&, const halfint&);
  friend const double operator/(const halfint&, const halfint&);
  friend const halfint operator*(int, const halfint&);
  friend const halfint operator*(const halfint&, int);
  friend const double operator*(double, const halfint&);
  friend const dcomplex operator*(dcomplex, const halfint&);
  friend const bool operator==(const halfint&, const halfint&);
  friend const bool operator!=(const halfint&, const halfint&);
  friend const bool operator<=(const halfint&, const halfint&);
  friend const bool operator>=(const halfint&, const halfint&);
  friend const bool operator<(const halfint&, const halfint&);
  friend const bool operator>(const halfint&, const halfint&);
  friend const int twice(const halfint&);
  operator double() const;
  operator dcomplex() const;
  friend const halfint habs(const halfint&);
  friend ostream& operator<<(ostream&, const halfint&);
  friend istream& operator>>(istream&, halfint&);
  friend bool whole(halfint);
  static const halfint _half;
  static const halfint _0;
  static const halfint _1;
  static const halfint _2;
  static const halfint _3;
  static const halfint _4;
  static const halfint _5;
  static const halfint _6;
  static const halfint _7;
  static const halfint _8;
  static const halfint _9;
  static const halfint _10;
private:
  /**
     Create a halfint by specifying twice its value.
     The second argument is a dummy, it's role is to differentiate from the constructor halfint(int),
     like in the case of i++ and ++i.
  */
  halfint(int _2x,int);
  int _2x;
};

extern const halfint& half;
extern const halfint& _0;
extern const halfint& _1;
extern const halfint& _2;
extern const halfint& _3;
extern const halfint& _4;
extern const halfint& _5;
extern const halfint& _6;
extern const halfint& _7;
extern const halfint& _8;
extern const halfint& _9;
extern const halfint& _10;

bool is_s1h(halfint h);
bool is_s3h(halfint h);
bool is_s5h(halfint h);
bool is_spin(halfint spin, halfint h);

#endif // HALFINT_HPP
