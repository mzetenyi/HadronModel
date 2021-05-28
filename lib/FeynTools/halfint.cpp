#include "halfint.hpp"
//#include <iostream>
//using std::cerr;
//using std::endl;

/** Default constructor, initialize to 0. */
halfint::halfint() : _2x(0) {}
/** Convert an int to halfint. The value of the resulting halfint will be x. */
halfint::halfint(int x) : _2x(2*x) {}
/**
   Create a halfint by specifying twice its value.
   The second argument is a dummy, it's role is to differentiate from the previous constructor,
   like in the case of i++ and ++i.
*/
halfint::halfint(int _2x,int) : _2x(_2x) {}
  /** copy constructor */
halfint::halfint(const halfint& h) : _2x(h._2x) {}

halfint& halfint::operator=(const halfint& h) {
  if (&h != this) {
    _2x = h._2x;
  }
  return *this;
}

halfint& halfint::operator=(int i) {
  _2x = 2*i;
  return *this;
}



const halfint& halfint::operator+() const { return *this; }
const halfint halfint::operator-() const { return halfint(-_2x,1); }
const halfint& halfint::operator++() {
  _2x += 2;
  return *this;
}
const halfint halfint::operator++(int) {
  halfint before(_2x,1);
  _2x += 2;
  return before;
}
const halfint& halfint::operator--() {
  _2x -= 2;
  return *this;
}
const halfint halfint::operator--(int) {
  halfint before(_2x,1);
  _2x -= 2;
  return before;
}
const halfint& halfint::operator+=(const halfint& h) {
  _2x += h._2x;
  return *this;
}
const halfint& halfint::operator-=(const halfint& h) {
  _2x -= h._2x;
  return *this;
}
const halfint& halfint::operator*=(int n) {
  _2x *= n;
  return *this;
}
const double halfint::operator*(double a) {
  return _2x/2. * a;
}
const double halfint::operator/(double a) {
  return _2x/2. / a;
}
const dcomplex halfint::operator*(dcomplex a) {
  return _2x/2. * a;
}
const dcomplex halfint::operator/(dcomplex a) {
  return _2x/2. / a;
}
const halfint operator+(const halfint& a, const halfint& b) {
  return halfint(a._2x + b._2x,1);
}
const halfint operator+(const halfint& a, int b) {
  return halfint(a._2x + 2*b,1);
}
const halfint operator+(int a, const halfint& b) {
  return halfint(2*a + b._2x,1);
}
double operator+(const halfint& a, double b) {
  return a._2x/2. + b;
}
double operator+(double a, const halfint& b) {
  return a + b._2x/2.;
}
dcomplex operator+(const halfint& a, dcomplex b) {
  return a._2x/2. + b;
}
dcomplex operator+(dcomplex a, const halfint& b) {
  return a + b._2x/2.;
}
const halfint operator-(const halfint& a,const halfint& b) {
  return halfint(a._2x - b._2x,1);
}
const halfint operator-(const halfint& a, int b) {
  return halfint(a._2x - 2*b,1);
}
const halfint operator-(int a, const halfint& b) {
  return halfint(2*a - b._2x,1);
}
double operator-(const halfint& a, double b) {
  return a._2x/2. - b;
}
double operator-(double a, const halfint& b) {
  return a - b._2x/2.;
}
dcomplex operator-(const halfint& a, dcomplex b) {
  return a._2x/2. - b;
}
dcomplex operator-(dcomplex a, const halfint& b) {
  return a - b._2x/2.;
}
const double operator*(const halfint& a,const halfint& b) {
  return (a._2x * b._2x)/4.;
}
const double operator/(const halfint& a,const halfint& b) {
  return a._2x/b._2x;
}
const halfint operator*(int n, const halfint& h) {
  return halfint(n*h._2x,1);
}
const halfint operator*(const halfint& h, int n) {
  return halfint(n*h._2x,1);
}
const double operator*(double a, const halfint& h){
  return h._2x/2. * a;
}
const dcomplex operator*(dcomplex a, const halfint& h){
  return h._2x/2. * a;
}
const bool operator==(const halfint& a, const halfint& b) {
  return a._2x == b._2x;
}
const bool operator!=(const halfint& a, const halfint& b) {
  return a._2x != b._2x;
}
const bool operator<=(const halfint& a, const halfint& b) {
  return a._2x <= b._2x;
}
const bool operator>=(const halfint& a, const halfint& b) {
  return a._2x >= b._2x;
}
const bool operator<(const halfint& a, const halfint& b) {
  return a._2x < b._2x;
}
const bool operator>(const halfint& a, const halfint& b) {
  return a._2x > b._2x;
}
const int twice(const halfint& h) {
  return h._2x;
}
halfint::operator double() const {
  return _2x/2.;
}
halfint::operator dcomplex() const {
  return dcomplex(_2x/2.);
}
const halfint habs(const halfint& h) {
  return (h>0 ? h : -h);
}

ostream& operator<<(ostream& out, const halfint& h) {
  if (h._2x % 2 == 0) {
    out << h._2x/2;
  } else {
    out << h._2x << "/2";
  }
  return out;
}

/**
   This is a very primitive implementation of operator>>, accepting only
   the form 'x/2' (no whitespaces etc.).
*/
istream& operator >> (istream& in, halfint& h) {
  h = halfint(0);
  halfint x(0);
  uint denom;
  char c;
  try {
    in >> x._2x;
    in.get(c);
    if (c!='/') {
      h = 2*x;
      return in;
    }
    in >> denom;
    if (denom!=2) { return in; }
    h = x;
  } catch(...) {}
  return in;
}

bool whole(halfint h) {
  return (h._2x%2 == 0);
}

const halfint halfint::_half(1,0);
const halfint halfint::_0(0,0);
const halfint halfint::_1(2,0);
const halfint halfint::_2(4,0);
const halfint halfint::_3(6,0);
const halfint halfint::_4(8,0);
const halfint halfint::_5(10,0);
const halfint halfint::_6(12,0);
const halfint halfint::_7(14,0);
const halfint halfint::_8(16,0);
const halfint halfint::_9(18,0);
const halfint halfint::_10(20,0);

const halfint& half(halfint::_half);
const halfint& _0(halfint::_0);
const halfint& _1(halfint::_1);
const halfint& _2(halfint::_2);
const halfint& _3(halfint::_3);
const halfint& _4(halfint::_4);
const halfint& _5(halfint::_5);
const halfint& _6(halfint::_6);
const halfint& _7(halfint::_7);
const halfint& _8(halfint::_8);
const halfint& _9(halfint::_9);
const halfint& _10(halfint::_10);

bool is_s1h(halfint h) {
  return ((h==-half) or (h==half));
}

bool is_s3h(halfint h) {
  return ((!whole(h)) and (-3*half<=h) and (h<=3*half));
}

bool is_s5h(halfint h) {
  return ((!whole(h)) and (-3*half<=h) and (h<=3*half));
}

bool is_spin(halfint spin, halfint h) {
  return ((whole(spin)==whole(h)) and (-spin<=h) and (h<=spin));
}
