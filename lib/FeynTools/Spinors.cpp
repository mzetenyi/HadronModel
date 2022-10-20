#include "Spinors.hpp"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

using namespace std;
using namespace Vectors;

namespace Spinors {

/// Create a nullspinor.
PauliSpinor::PauliSpinor() : x1(0), x2(0) {}
/// Copy constructor.
PauliSpinor::PauliSpinor(const PauliSpinor& s) : x1(s.x1), x2(s.x2) {}
/// Create a spinor with the given components.
PauliSpinor::PauliSpinor(dcomplex x1, dcomplex x2) : x1(x1), x2(x2) {}
/// Destructor. Does nothing.
PauliSpinor::~PauliSpinor() {}
/// Assignment operator.
PauliSpinor& PauliSpinor::operator=(const PauliSpinor& s) {
  x1 = s.x1;
  x2 = s.x2;
  return *this;
}
/// Adjoint
AdPauliSpinor PauliSpinor::adj() const {
  return AdPauliSpinor(conj(x1), conj(x2));
}
AdPauliSpinor adj(const PauliSpinor& s) {
  return AdPauliSpinor(conj(s.x1), conj(s.x2));
}
/// Positive sign.
const PauliSpinor& PauliSpinor::operator+() const { return *this; }
/// Negative sign.
PauliSpinor PauliSpinor::operator-() const { return PauliSpinor(-x1, -x2); }
/// Spinor addition.
PauliSpinor PauliSpinor::operator+(const PauliSpinor& s) const {
  return PauliSpinor(x1 + s.x1, x2 + s.x2);
}
/// Spinor subtraction.
PauliSpinor PauliSpinor::operator-(const PauliSpinor& s) const {
  return PauliSpinor(x1 - s.x1, x2 - s.x2);
}
/// Increment operator.
PauliSpinor& PauliSpinor::operator+=(const PauliSpinor& s) {
  x1 += s.x1;
  x2 += s.x2;
  return *this;
}
/// Decrement operator.
PauliSpinor& PauliSpinor::operator-=(const PauliSpinor& s) {
  x1 -= s.x1;
  x2 -= s.x2;
  return *this;
}
/// Multiplication with scalar (dcomplex).
PauliSpinor PauliSpinor::operator*(const dcomplex& z) const {
  return PauliSpinor(z * x1, z * x2);
}
/// Division by scalar.
PauliSpinor PauliSpinor::operator/(const dcomplex& z) const {
  return PauliSpinor(x1 / z, x2 / z);
}

PauliSpinor& PauliSpinor::operator*=(const dcomplex& z) {
  x1 *= z;
  x2 *= z;
  return *this;
}

PauliSpinor& PauliSpinor::operator/=(const dcomplex& z) {
  x1 /= z;
  x2 /= z;
  return *this;
}
/// Friend operator for scalar*spinor.
PauliSpinor operator*(const dcomplex& z, const PauliSpinor& s) {
  return PauliSpinor(z * s.x1, z * s.x2);
}
/// Matrix*spinor product.
PauliSpinor operator*(const PauliMatrix& tau, const PauliSpinor& s) {
  return PauliSpinor(tau.x11 * s.x1 + tau.x12 * s.x2,
                     tau.x21 * s.x1 + tau.x22 * s.x2);
}
/// Scalar product.
dcomplex operator*(const AdPauliSpinor& a, const PauliSpinor& s) {
  return a.x1 * s.x1 + a.x2 * s.x2;
}
/// Diadic product.
PauliMatrix operator*(const PauliSpinor& s, const AdPauliSpinor& a) {
  return PauliMatrix(s.x1 * a.x1, s.x1 * a.x2, s.x2 * a.x1, s.x2 * a.x2);
}
/// Equality check.
bool PauliSpinor::operator==(const PauliSpinor& s) const {
  return (x1 == s.x1 and x2 == s.x2);
}
/// Unequality check.
bool PauliSpinor::operator!=(const PauliSpinor& s) const {
  return (x1 != s.x1 or x2 != s.x2);
}
/// Print the spinor to the given stream in the form '(x,y)'.
ostream& operator<<(ostream& out, const PauliSpinor& s) {
  out << "(" << s.x1 << ", " << s.x2 << ")";
  return out;
}
/**
   Read the spinor from the given stream. It must be in the form '[x,y]'.
   This is a very primitive implementation of operator>>, accepting only
   the form '[vx,vy]' (no whitespaces etc.). [] are used because () is
   interpreted by bash. Must be improved. (A separate parser class should
   be used.)

   \todo This is a very primitive implementation of operator>>, accepting only
   the form '[vx,vy]' (no whitespaces etc.). [] are used because () is
   interpreted by bash. Must be improved. (A separate parser class should
   be used.)
*/
istream& operator>>(istream& in, PauliSpinor& s) {
  s = PauliSpinor(0, 0);
  PauliSpinor ss(0, 0);
  char c;
  try {
    in.get(c);
    if (c != '[') {
      return in;
    }
    in >> ss.x1;
    in.get(c);
    if (c != ',') {
      return in;
    }
    in >> ss.x2;
    in.get(c);
    if (c != ']') {
      return in;
    }
    s = ss;  // put result in s ONLY if everything was OK
  } catch (...) {
  }
  return in;
}

/// Decide if any component is NaN.
bool PauliSpinor::isnan() const { return ::isnan(x1) or ::isnan(x2); }

const PauliSpinor Khi_plus(1, 0);
const PauliSpinor Khi_minus(0, 1);
const PauliSpinor Khi_null(0, 0);

/**
   Return the Pauli spinor of the given spin.
   (The value of the argument is the spin.)
*/
const PauliSpinor& Khi(halfint spin) {
  assert(spin == -half or spin == half);
  if (spin == -half) {
    return Khi_minus;
  }
  if (spin == half) {
    return Khi_plus;
  }
  return Khi_null;  // never get here
}

/// Create a nullspinor.
AdPauliSpinor::AdPauliSpinor() : x1(0), x2(0) {}
/// Copy constructor.
AdPauliSpinor::AdPauliSpinor(const AdPauliSpinor& a) : x1(a.x1), x2(a.x2) {}
/// Create a spinor with the given components.
AdPauliSpinor::AdPauliSpinor(dcomplex x1, dcomplex x2) : x1(x2), x2(x2) {}
/// Destructor. Does nothing.
AdPauliSpinor::~AdPauliSpinor() {}
/// Assignment operator.
AdPauliSpinor& AdPauliSpinor::operator=(const AdPauliSpinor& a) {
  x1 = a.x1;
  x2 = a.x2;
  return *this;
}
/// Adjoint
PauliSpinor AdPauliSpinor::adj() const {
  return PauliSpinor(conj(x1), conj(x2));
}
PauliSpinor adj(const AdPauliSpinor& a) {
  return PauliSpinor(conj(a.x1), conj(a.x2));
}
/// Positive sign.
const AdPauliSpinor& AdPauliSpinor::operator+() const { return *this; }
/// Negative sign.
AdPauliSpinor AdPauliSpinor::operator-() const {
  return AdPauliSpinor(-x1, -x2);
}
/// Spinor addition.
AdPauliSpinor AdPauliSpinor::operator+(const AdPauliSpinor& a) const {
  return AdPauliSpinor(x1 + a.x1, x2 + a.x2);
}
/// Spinor subtraction.
AdPauliSpinor AdPauliSpinor::operator-(const AdPauliSpinor& a) const {
  return AdPauliSpinor(x1 - a.x1, x2 - a.x2);
}
/// Increment operator.
AdPauliSpinor& AdPauliSpinor::operator+=(const AdPauliSpinor& a) {
  x1 += a.x1;
  x2 += a.x2;
  return *this;
}
/// Decrement operator.
AdPauliSpinor& AdPauliSpinor::operator-=(const AdPauliSpinor& a) {
  x1 -= a.x1;
  x2 -= a.x2;
  return *this;
}
/// Multiplication with scalar (dcomplex).
AdPauliSpinor AdPauliSpinor::operator*(const dcomplex& z) const {
  return AdPauliSpinor(z * x1, z * x2);
}
/// Division by scalar.
AdPauliSpinor AdPauliSpinor::operator/(const dcomplex& z) const {
  return AdPauliSpinor(x1 / z, x2 / z);
}

AdPauliSpinor& AdPauliSpinor::operator*=(const dcomplex& z) {
  x1 *= z;
  x2 *= z;
  return *this;
}

AdPauliSpinor& AdPauliSpinor::operator/=(const dcomplex& z) {
  x1 /= z;
  x2 /= z;
  return *this;
}
/// Friend operator for scalar*spinor.
AdPauliSpinor operator*(const dcomplex& z, const AdPauliSpinor& a) {
  return AdPauliSpinor(z * a.x1, z * a.x2);
}
/// Spinor*matrix product.
AdPauliSpinor operator*(const AdPauliSpinor& a, const PauliMatrix& tau) {
  return AdPauliSpinor(a.x1 * tau.x11 + a.x2 * tau.x21,
                       a.x1 * tau.x12 + a.x2 * tau.x22);
}
/// Equality check.
bool AdPauliSpinor::operator==(const AdPauliSpinor& a) const {
  return (x1 == a.x1 and x2 == a.x2);
}
/// Unequality check.
bool AdPauliSpinor::operator!=(const AdPauliSpinor& a) const {
  return (x1 != a.x1 or x2 != a.x2);
}
/// Print the spinor to the given stream in the form '(x,y)'.
ostream& operator<<(ostream& out, const AdPauliSpinor& a) {
  out << "(" << a.x1 << ", " << a.x2 << ")";
  return out;
}
/**
   Read the spinor from the given stream. It must be in the form '[x,y]'.
   This is a very primitive implementation of operator>>, accepting only
   the form '[vx,vy]' (no whitespaces etc.). [] are used because () is
   interpreted by bash. Must be improved. (A separate parser class should
   be used.)

   \todo This is a very primitive implementation of operator>>, accepting only
   the form '[vx,vy]' (no whitespaces etc.). [] are used because () is
   interpreted by bash. Must be improved. (A separate parser class should
   be used.)
*/
istream& operator>>(istream& in, AdPauliSpinor& a) {
  a = AdPauliSpinor(0, 0);
  AdPauliSpinor aa(0, 0);
  char c;
  try {
    in.get(c);
    if (c != '[') {
      return in;
    }
    in >> aa.x1;
    in.get(c);
    if (c != ',') {
      return in;
    }
    in >> aa.x2;
    in.get(c);
    if (c != ']') {
      return in;
    }
    a = aa;  // put result in s ONLY if everything was OK
  } catch (...) {
  }
  return in;
}

/// Decide if any component is NaN.
bool AdPauliSpinor::isnan() const { return ::isnan(x1) or ::isnan(x2); }

/// Create a nullmatrix.
PauliMatrix::PauliMatrix() : x11(0), x12(0), x21(0), x22(0) {}
/// Copy constructor.
PauliMatrix::PauliMatrix(const PauliMatrix& _t)
    : x11(_t.x11), x12(_t.x12), x21(_t.x21), x22(_t.x22) {}
/// Create a Pauli matrix from components.
PauliMatrix::PauliMatrix(dcomplex x11, dcomplex x12, dcomplex x21, dcomplex x22)
    : x11(x11), x12(x12), x21(x21), x22(x22) {
  /*
  cerr << "PauliMatrix::constr" << endl;
  PR(x11);PR(x12);PR(x21);PR(x22);
  PR(*this);
  */
}
/// Destructor. Does nothing.
PauliMatrix::~PauliMatrix() {}
/// Assignment operator.
PauliMatrix& PauliMatrix::operator=(const PauliMatrix& _t) {
  x11 = _t.x11;
  x12 = _t.x12;
  x21 = _t.x21;
  x22 = _t.x22;
  return *this;
}
/// Adjoint
PauliMatrix PauliMatrix::adj() const {
  return PauliMatrix(std::conj(x11), std::conj(x21), std::conj(x12),
                     std::conj(x22));
}
PauliMatrix adj(const PauliMatrix& _t) {
  return PauliMatrix(std::conj(_t.x11), std::conj(_t.x21), std::conj(_t.x12),
                     std::conj(_t.x22));
}
/// Positive sign.
const PauliMatrix& PauliMatrix::operator+() const { return *this; }
/// Negative sign.
PauliMatrix PauliMatrix::operator-() const {
  return PauliMatrix(-x11, -x12, -x21, -x22);
}
/// Matrix addition.
PauliMatrix PauliMatrix::operator+(const PauliMatrix& _t) const {
  return PauliMatrix(x11 + _t.x11, x12 + _t.x12, x21 + _t.x21, x22 + _t.x22);
}
/// Matrix subtraction.
PauliMatrix PauliMatrix::operator-(const PauliMatrix& _t) const {
  return PauliMatrix(x11 - _t.x11, x12 - _t.x12, x21 - _t.x21, x22 - _t.x22);
}
/// Increment operator.
PauliMatrix& PauliMatrix::operator+=(const PauliMatrix& _t) {
  x11 += _t.x11;
  x12 += _t.x12;
  x21 += _t.x21;
  x22 += _t.x22;
  return *this;
}
/// Decrement operator.
PauliMatrix& PauliMatrix::operator-=(const PauliMatrix& _t) {
  x11 -= _t.x11;
  x12 -= _t.x12;
  x21 -= _t.x21;
  x22 -= _t.x22;
  return *this;
}
/// Matrix product.
PauliMatrix operator*(const PauliMatrix& _t1, const PauliMatrix& _t2) {
  return PauliMatrix(_t1.x11 * _t2.x11 + _t1.x12 * _t2.x21,
                     _t1.x11 * _t2.x12 + _t1.x12 * _t2.x22,
                     _t1.x21 * _t2.x11 + _t1.x22 * _t2.x21,
                     _t1.x21 * _t2.x12 + _t1.x22 * _t2.x22);
}
/// Multiplication with scalar (dcomplex).
PauliMatrix PauliMatrix::operator*(const dcomplex& z) const {
  return PauliMatrix(x11 * z, x12 * z, x21 * z, x22 * z);
}
/// Division by scalar.
PauliMatrix PauliMatrix::operator/(dcomplex z) const {
  // cerr << "PauliMatrix::operator / (dcomplex z) " << z << endl;
  return PauliMatrix(x11 / z, x12 / z, x21 / z, x22 / z);
}
PauliMatrix& PauliMatrix::operator*=(const dcomplex& z) {
  x11 *= z;
  x12 *= z;
  x21 *= z;
  x22 *= z;
  return *this;
}
PauliMatrix& PauliMatrix::operator/=(const dcomplex& z) {
  x11 /= z;
  x12 /= z;
  x21 /= z;
  x22 /= z;
  return *this;
}
/// Friend operator for scalar*matrix.
PauliMatrix operator*(const dcomplex& z, const PauliMatrix& _t) {
  return PauliMatrix(z * _t.x11, z * _t.x12, z * _t.x21, z * _t.x22);
}
/// Equality check.
bool PauliMatrix::operator==(const PauliMatrix& _t) const {
  return x11 == _t.x11 and x12 == _t.x12 and x21 == _t.x21 and x22 == _t.x22;
}
/// Unequality check.
bool PauliMatrix::operator!=(const PauliMatrix& _t) const {
  return x11 != _t.x11 or x12 != _t.x12 or x21 != _t.x21 or x22 != _t.x22;
}
/// Print the matrix to the given stream in the form '((x11,x12),(...))'.
ostream& operator<<(ostream& out, const PauliMatrix& _t) {
  out << "((" << _t.x11 << "," << _t.x12 << "), (" << _t.x21 << "," << _t.x22
      << "))";
  return out;
}
/// Print to a stream in the form of a table.
void PauliMatrix::output(ostream& out, int fieldWidth) const {
  out << "( " << setw(fieldWidth) << x11 << "  " << setw(fieldWidth) << x12
      << " )" << endl;
  out << "( " << setw(fieldWidth) << x21 << "  " << setw(fieldWidth) << x22
      << " )" << endl;
}
/// Trace.
dcomplex PauliMatrix::trace() const { return x11 + x22; }

PauliMatrix PauliMatrix::conj() const {
  return PauliMatrix(std::conj(x11), std::conj(x12), std::conj(x21),
                     std::conj(x22));
}
PauliMatrix PauliMatrix::transpose() const {
  return PauliMatrix(x11, x21, x12, x22);
}

const PauliMatrix tau1(0, 1, 1, 0);
const PauliMatrix tau2(0, -i_, i_, 0);
const PauliMatrix tau3(1, 0, 0, -1);
const PauliMatrix tau_null(0, 0, 0, 0);
const PauliMatrix tau_unit(1, 0, 0, 1);

const Tau_i& Tau_i::instance() { return tau_; }

PauliMatrix Tau_i::operator()(int i) const {
  assert(1 <= i and i < 4);
  if (i == 1) {
    return tau1;
  }
  if (i == 2) {
    return tau2;
  }
  if (i == 3) {
    return tau3;
  }
  return tau_null;  // never get here
}

const PauliMatrix& Tau_i::operator[](int i) const {
  assert(1 <= i and i < 4);
  if (i == 1) {
    return tau1;
  }
  if (i == 2) {
    return tau2;
  }
  if (i == 3) {
    return tau3;
  }
  return tau_null;  // never get here
}

Tau_i::Tau_i() {}

const Tau_i& tau_ = Tau_i::instance();

const Tau_i Tau_i::tau_;

/// Create a nullspinor.
DiracSpinor::DiracSpinor() {
  for (int i(0); i < 4; i++) {
    x[i] = 0;
  }
}
/// Copy constructor.
DiracSpinor::DiracSpinor(const DiracSpinor& u) {
  for (int i(0); i < 4; i++) {
    x[i] = u.x[i];
  }
}
/// Create a spinor with the given components.
DiracSpinor::DiracSpinor(dcomplex x0, dcomplex x1, dcomplex x2, dcomplex x3) {
  x[0] = x0;
  x[1] = x1;
  x[2] = x2;
  x[3] = x3;
}
/// Create a spinor from two Pauli spinors.
DiracSpinor::DiracSpinor(const PauliSpinor& s1, const PauliSpinor& s2) {
  x[0] = s1.x1;
  x[1] = s1.x2;
  x[2] = s2.x1;
  x[3] = s2.x2;
}
/// Destructor. Does nothing.
DiracSpinor::~DiracSpinor() {}
/// Assignment operator.
DiracSpinor& DiracSpinor::operator=(const DiracSpinor& u) {
  x[0] = u.x[0];
  x[1] = u.x[1];
  x[2] = u.x[2];
  x[3] = u.x[3];
  return *this;
}
/// Positive sign.
const DiracSpinor& DiracSpinor::operator+() const { return *this; }
/// Negative sign.
DiracSpinor DiracSpinor::operator-() const { return -1 * (*this); }
/// Spinor addition.
DiracSpinor DiracSpinor::operator+(const DiracSpinor& u) const {
  return DiracSpinor(x[0] + u.x[0], x[1] + u.x[1], x[2] + u.x[2],
                     x[3] + u.x[3]);
}
/// Spinor subtraction.
DiracSpinor DiracSpinor::operator-(const DiracSpinor& u) const {
  return DiracSpinor(x[0] - u.x[0], x[1] - u.x[1], x[2] - u.x[2],
                     x[3] - u.x[3]);
}
/// Increment operator.
DiracSpinor& DiracSpinor::operator+=(const DiracSpinor& u) {
  for (int i(0); i < 4; i++) {
    x[i] += u.x[i];
  }
  return *this;
}
/// Decrement operator.
DiracSpinor& DiracSpinor::operator-=(const DiracSpinor& u) {
  for (int i(0); i < 4; i++) {
    x[i] -= u.x[i];
  }
  return *this;
}
/// Multiplication with scalar (dcomplex).
DiracSpinor DiracSpinor::operator*(const dcomplex& z) const {
  return DiracSpinor(x[0] * z, x[1] * z, x[2] * z, x[3] * z);
}
/// Division by scalar.
DiracSpinor DiracSpinor::operator/(const dcomplex& z) const {
  return DiracSpinor(x[0] / z, x[1] / z, x[2] / z, x[3] / z);
}
DiracSpinor& DiracSpinor::operator*=(const dcomplex& z) {
  for (int i(0); i < 4; i++) {
    x[i] *= z;
  }
  return *this;
}
DiracSpinor& DiracSpinor::operator/=(const dcomplex& z) {
  for (int i(0); i < 4; i++) {
    x[i] /= z;
  }
  return *this;
}
/// Friend operator for scalar*spinor.
DiracSpinor operator*(const dcomplex& z, const DiracSpinor& u) {
  return DiracSpinor(u.x[0] * z, u.x[1] * z, u.x[2] * z, u.x[3] * z);
}
/// Matrix*spinor product.
DiracSpinor operator*(const DiracMatrix& gamma, const DiracSpinor& u) {
  DiracSpinor v(0, 0, 0, 0);
  for (int i(0); i < 4; i++) {
    for (int j(0); j < 4; j++) {
      v.x[i] += gamma.x[i][j] * u.x[j];
    }
  }
  return v;
}
/// Scalar product.
dcomplex operator*(const AdDiracSpinor& ubar, const DiracSpinor& u) {
  dcomplex z(0);
  for (int i(0); i < 4; i++) {
    z += ubar.x[i] * u.x[i];
  }
  return z;
}
/// Diadic product.
DiracMatrix operator*(const DiracSpinor& u, const AdDiracSpinor& ubar) {
  DiracMatrix gamma;
  for (int i(0); i < 4; i++) {
    for (int j(0); j < 4; j++) {
      gamma.x[i][j] = u.x[i] * ubar.x[j];
    }
  }
  return gamma;
}
/// Equality check.
bool DiracSpinor::operator==(const DiracSpinor& u) const {
  for (int i(0); i < 4; i++) {
    if (x[i] != u.x[i]) {
      return false;
    }
  }
  return true;
}
/// Unequality check.
bool DiracSpinor::operator!=(const DiracSpinor& u) const {
  for (int i(0); i < 4; i++) {
    if (x[i] != u.x[i]) {
      return true;
    }
  }
  return false;
}
/// Print the spinor to the given stream in the form '(x1,x2,x3,x4)'.
ostream& operator<<(ostream& out, const DiracSpinor& u) {
  out << "(" << u.x[0] << ", " << u.x[1] << ", " << u.x[2] << ", " << u.x[3]
      << ")";
  return out;
}
/// Read the spinor from the given stream. It must be in the form
/// '[x1,x2,x3,x4]'.
istream& operator>>(istream& in, DiracSpinor& u) {
  u = DiracSpinor(0, 0, 0, 0);
  DiracSpinor v(0, 0, 0, 0);
  char c;
  try {
    in.get(c);
    if (c != '[') {
      return in;
    }
    in >> v.x[0];
    in.get(c);
    if (c != ',') {
      return in;
    }
    in >> v.x[1];
    in.get(c);
    if (c != ',') {
      return in;
    }
    in >> v.x[2];
    in.get(c);
    if (c != ',') {
      return in;
    }
    in >> v.x[3];
    in.get(c);
    if (c != ']') {
      return in;
    }
    u = v;  // put result in s ONLY if everything was OK
  } catch (...) {
  }
  return in;
}
/// Decide if any component is NaN.
bool DiracSpinor::isnan() const {
  for (int i(0); i < 4; i++) {
    if (::isnan(x[i])) {
      return true;
    }
  }
  return false;
}

AdDiracSpinor DiracSpinor::adj() const {
  return AdDiracSpinor(conj(x[0]), conj(x[1]), conj(x[2]), conj(x[3])) *
         gamma0_;
}

AdDiracSpinor adj(const DiracSpinor& u) {
  return AdDiracSpinor(conj(u.x[0]), conj(u.x[1]), conj(u.x[2]), conj(u.x[3])) *
         gamma0_;
}

DiracSpinor U_(const FourVector& p, halfint spin) {
  assert(spin == -half or spin == half);
  // cerr << "DiracSpinor u_(const FourVector& p, int spin)" << endl;
  // PR(p); PR(p*p);
  double m2 = p * p;
  if (m2 < 0) {
    if (-m2 / (p(0) * p(0)) > 0.000001) {
      cerr << "spacelike fourmomentum in u_: " << p << endl;
      exit(0);
    }
    m2 = 0;
  }
  double m = sqrt(m2);
  PauliMatrix tau_p(tau_null);
  for (int i(1); i < 4; i++) {
    tau_p += tau_(i) * p(i);
  }
  double norm = sqrt(p(0) + m);
  // PR(norm);
  PauliSpinor khi = Khi(spin);
  // PR(khi);
  return DiracSpinor(norm * khi, 1. / norm * (tau_p * khi));
}

DiracSpinor V_(const FourVector& p, halfint spin) {
  assert(spin == -half or spin == half);
  double m2 = p * p;
  if (m2 < 0) {
    if (-m2 / (p(0) * p(0)) > 0.000001) {
      cerr << "spacelike fourmomentum in u_: " << p << endl;
      exit(0);
    }
    m2 = 0;
  }
  double m = sqrt(m2);
  PauliMatrix tau_p(tau_null);
  for (int i(1); i < 4; i++) {
    tau_p += tau_(i) * p(i);
  }
  double norm = sqrt(p(0) + m);
  PauliSpinor khi = Khi(-spin);
  return DiracSpinor(spin / norm * (tau_p * khi), spin * norm * khi);
}

/// Create a nullspinor.
AdDiracSpinor::AdDiracSpinor() {
  for (int i(0); i < 4; i++) {
    x[i] = 0;
  }
}
/// Copy constructor.
AdDiracSpinor::AdDiracSpinor(const AdDiracSpinor& u) {
  for (int i(0); i < 4; i++) {
    x[i] = u.x[i];
  }
}
/// Create a spinor with the given components.
AdDiracSpinor::AdDiracSpinor(dcomplex x0, dcomplex x1, dcomplex x2,
                             dcomplex x3) {
  x[0] = x0;
  x[1] = x1;
  x[2] = x2;
  x[3] = x3;
}
/// Create a spinor from two Pauli spinors.
AdDiracSpinor::AdDiracSpinor(const AdPauliSpinor& s1, const AdPauliSpinor& s2) {
  x[0] = s1.x1;
  x[1] = s1.x2;
  x[2] = s2.x1;
  x[3] = s2.x2;
}
/// Destructor. Does nothing.
AdDiracSpinor::~AdDiracSpinor() {}
/// Assignment operator.
AdDiracSpinor& AdDiracSpinor::operator=(const AdDiracSpinor& u) {
  x[0] = u.x[0];
  x[1] = u.x[1];
  x[2] = u.x[2];
  x[3] = u.x[3];
  return *this;
}
/// Positive sign.
const AdDiracSpinor& AdDiracSpinor::operator+() const { return *this; }
/// Negative sign.
AdDiracSpinor AdDiracSpinor::operator-() const { return -1 * (*this); }
/// Spinor addition.
AdDiracSpinor AdDiracSpinor::operator+(const AdDiracSpinor& u) const {
  return AdDiracSpinor(x[0] + u.x[0], x[1] + u.x[1], x[2] + u.x[2],
                       x[3] + u.x[3]);
}
/// Spinor subtraction.
AdDiracSpinor AdDiracSpinor::operator-(const AdDiracSpinor& u) const {
  return AdDiracSpinor(x[0] - u.x[0], x[1] - u.x[1], x[2] - u.x[2],
                       x[3] - u.x[3]);
}
/// Increment operator.
AdDiracSpinor& AdDiracSpinor::operator+=(const AdDiracSpinor& u) {
  for (int i(0); i < 4; i++) {
    x[i] += u.x[i];
  }
  return *this;
}
/// Decrement operator.
AdDiracSpinor& AdDiracSpinor::operator-=(const AdDiracSpinor& u) {
  for (int i(0); i < 4; i++) {
    x[i] -= u.x[i];
  }
  return *this;
}
/// Multiplication with scalar (dcomplex).
AdDiracSpinor AdDiracSpinor::operator*(const dcomplex& z) const {
  return AdDiracSpinor(x[0] * z, x[1] * z, x[2] * z, x[3] * z);
}
/// Division by scalar.
AdDiracSpinor AdDiracSpinor::operator/(const dcomplex& z) const {
  return AdDiracSpinor(x[0] / z, x[1] / z, x[2] / z, x[3] / z);
}
AdDiracSpinor& AdDiracSpinor::operator*=(const dcomplex& z) {
  for (int i(0); i < 4; i++) {
    x[i] *= z;
  }
  return *this;
}
AdDiracSpinor& AdDiracSpinor::operator/=(const dcomplex& z) {
  for (int i(0); i < 4; i++) {
    x[i] /= z;
  }
  return *this;
}
/// Friend operator for scalar*spinor.
AdDiracSpinor operator*(const dcomplex& z, const AdDiracSpinor& u) {
  return AdDiracSpinor(u.x[0] * z, u.x[1] * z, u.x[2] * z, u.x[3] * z);
}
/// spinor*matrix product.
AdDiracSpinor operator*(const AdDiracSpinor& u, const DiracMatrix& gamma) {
  // cerr << "AdDiracSpinor operator * (const AdDiracSpinor& u, const
  // DiracMatrix& gamma)" << endl;
  AdDiracSpinor v(0, 0, 0, 0);
  for (int i(0); i < 4; i++) {
    for (int j(0); j < 4; j++) {
      //        PR(i); PR(j); PR(u.x[j]); PR(gamma.x[j][i]);
      v.x[i] += u.x[j] * gamma.x[j][i];
    }
  }
  return v;
}
/// Equality check.
bool AdDiracSpinor::operator==(const AdDiracSpinor& u) const {
  for (int i(0); i < 4; i++) {
    if (x[i] != u.x[i]) {
      return false;
    }
  }
  return true;
}
/// Unequality check.
bool AdDiracSpinor::operator!=(const AdDiracSpinor& u) const {
  for (int i(0); i < 4; i++) {
    if (x[i] != u.x[i]) {
      return true;
    }
  }
  return false;
}
/// Print the spinor to the given stream in the form '(x1,x2,x3,x4)'.
ostream& operator<<(ostream& out, const AdDiracSpinor& u) {
  out << "(" << u.x[0] << ", " << u.x[1] << ", " << u.x[2] << ", " << u.x[3]
      << ")";
  return out;
}
/// Read the spinor from the given stream. It must be in the form
/// '[x1,x2,x3,x4]'.
istream& operator>>(istream& in, AdDiracSpinor& u) {
  u = AdDiracSpinor(0, 0, 0, 0);
  AdDiracSpinor v(0, 0, 0, 0);
  char c;
  try {
    in.get(c);
    if (c != '[') {
      return in;
    }
    in >> v.x[0];
    in.get(c);
    if (c != ',') {
      return in;
    }
    in >> v.x[1];
    in.get(c);
    if (c != ',') {
      return in;
    }
    in >> v.x[2];
    in.get(c);
    if (c != ',') {
      return in;
    }
    in >> v.x[3];
    in.get(c);
    if (c != ']') {
      return in;
    }
    u = v;  // put result in s ONLY if everything was OK
  } catch (...) {
  }
  return in;
}
/// Decide if any component is NaN.
bool AdDiracSpinor::isnan() const {
  for (int i(0); i < 4; i++) {
    if (::isnan(x[i])) {
      return true;
    }
  }
  return false;
}

DiracSpinor AdDiracSpinor::adj() const {
  return gamma0_ * DiracSpinor(conj(x[0]), conj(x[1]), conj(x[2]), conj(x[3]));
}

DiracSpinor adj(const AdDiracSpinor& ubar) {
  return gamma0_ * DiracSpinor(conj(ubar.x[0]), conj(ubar.x[1]),
                               conj(ubar.x[2]), conj(ubar.x[3]));
}

AdDiracSpinor Ubar_(const FourVector& p, halfint spin) {
  return adj(U_(p, spin));
}
AdDiracSpinor Vbar_(const FourVector& p, halfint spin) {
  return adj(V_(p, spin));
}
/// Create a nullmatrix.
DiracMatrix::DiracMatrix() {
  for (int i(0); i < 4; i++) {
    for (int j(0); j < 4; j++) {
      x[i][j] = 0;
    }
  }
}
/// Copy constructor.
DiracMatrix::DiracMatrix(const DiracMatrix& g) {
  for (int i(0); i < 4; i++) {
    for (int j(0); j < 4; j++) {
      x[i][j] = g.x[i][j];
    }
  }
}
/// Create a Dirac matrix from components.
DiracMatrix::DiracMatrix(dcomplex g00, dcomplex g01, dcomplex g02, dcomplex g03,
                         dcomplex g10, dcomplex g11, dcomplex g12, dcomplex g13,
                         dcomplex g20, dcomplex g21, dcomplex g22, dcomplex g23,
                         dcomplex g30, dcomplex g31, dcomplex g32,
                         dcomplex g33) {
  x[0][0] = g00;
  x[0][1] = g01;
  x[0][2] = g02;
  x[0][3] = g03;
  x[1][0] = g10;
  x[1][1] = g11;
  x[1][2] = g12;
  x[1][3] = g13;
  x[2][0] = g20;
  x[2][1] = g21;
  x[2][2] = g22;
  x[2][3] = g23;
  x[3][0] = g30;
  x[3][1] = g31;
  x[3][2] = g32;
  x[3][3] = g33;
}

/// Create a Dirac matrix from Pauli matrices.
DiracMatrix::DiracMatrix(PauliMatrix t1, PauliMatrix t2, PauliMatrix t3,
                         PauliMatrix t4) {
  x[0][0] = t1.x11;
  x[0][1] = t1.x12;
  x[0][2] = t2.x11;
  x[0][3] = t2.x12;
  x[1][0] = t1.x21;
  x[1][1] = t1.x22;
  x[1][2] = t2.x21;
  x[1][3] = t2.x22;
  x[2][0] = t3.x11;
  x[2][1] = t3.x12;
  x[2][2] = t4.x11;
  x[2][3] = t4.x12;
  x[3][0] = t3.x21;
  x[3][1] = t3.x22;
  x[3][2] = t4.x21;
  x[3][3] = t4.x22;
}
/// Destructor. Does nothing.
DiracMatrix::~DiracMatrix() {}
/// Assignment operator.
DiracMatrix& DiracMatrix::operator=(const DiracMatrix& g) {
  for (int i(0); i < 4; i++) {
    for (int j(0); j < 4; j++) {
      x[i][j] = g.x[i][j];
    }
  }
  return *this;
}
/// Adjoint
DiracMatrix DiracMatrix::adj() const {
  return gamma0_ * (*this).trans_conj() * gamma0_;
}
DiracMatrix adj(const DiracMatrix& g) {
  return gamma0_ * g.trans_conj() * gamma0_;
}
/// Positive sign.
const DiracMatrix& DiracMatrix::operator+() const { return *this; }
/// Negative sign.
DiracMatrix DiracMatrix::operator-() const { return (-1) * (*this); }
/// Matrix addition.
DiracMatrix DiracMatrix::operator+(const DiracMatrix& g) const {
  DiracMatrix gg;
  for (int i(0); i < 4; i++) {
    for (int j(0); j < 4; j++) {
      gg.x[i][j] = x[i][j] + g.x[i][j];
    }
  }
  return gg;
}
/// Matrix subtraction.
DiracMatrix DiracMatrix::operator-(const DiracMatrix& g) const {
  DiracMatrix gg;
  for (int i(0); i < 4; i++) {
    for (int j(0); j < 4; j++) {
      gg.x[i][j] = x[i][j] - g.x[i][j];
    }
  }
  return gg;
}
/// Increment DiracMatrix::operator.
DiracMatrix& DiracMatrix::operator+=(const DiracMatrix& g) {
  for (int i(0); i < 4; i++) {
    for (int j(0); j < 4; j++) {
      x[i][j] += g.x[i][j];
    }
  }
  return *this;
}
/// Decrement operator.
DiracMatrix& DiracMatrix::operator-=(const DiracMatrix& g) {
  for (int i(0); i < 4; i++) {
    for (int j(0); j < 4; j++) {
      x[i][j] -= g.x[i][j];
    }
  }
  return *this;
}
/// Matrix product.
DiracMatrix operator*(const DiracMatrix& g1, const DiracMatrix& g2) {
  DiracMatrix g;
  for (int i(0); i < 4; i++) {
    for (int j(0); j < 4; j++) {
      for (int k(0); k < 4; k++) {
        g.x[i][j] += g1.x[i][k] * g2.x[k][j];
      }
    }
  }
  return g;
}
/// Multiplication with scalar (dcomplex).
DiracMatrix DiracMatrix::operator*(const dcomplex& z) const {
  DiracMatrix g;
  for (int i(0); i < 4; i++) {
    for (int j(0); j < 4; j++) {
      g.x[i][j] = x[i][j] * z;
    }
  }
  return g;
}
/// Division by scalar.
DiracMatrix DiracMatrix::operator/(const dcomplex& z) const {
  DiracMatrix g;
  for (int i(0); i < 4; i++) {
    for (int j(0); j < 4; j++) {
      g.x[i][j] = x[i][j] / z;
    }
  }
  return g;
}
DiracMatrix& DiracMatrix::operator*=(const dcomplex& z) {
  for (int i(0); i < 4; i++) {
    for (int j(0); j < 4; j++) {
      x[i][j] *= z;
    }
  }
  return *this;
}
DiracMatrix& DiracMatrix::operator/=(const dcomplex& z) {
  for (int i(0); i < 4; i++) {
    for (int j(0); j < 4; j++) {
      x[i][j] /= z;
    }
  }
  return *this;
}
/// Friend operator for scalar*matrix.
DiracMatrix operator*(const dcomplex& z, const DiracMatrix& g) {
  DiracMatrix gg;
  for (int i(0); i < 4; i++) {
    for (int j(0); j < 4; j++) {
      gg.x[i][j] = z * g.x[i][j];
    }
  }
  return gg;
}
/// Equality check.
bool DiracMatrix::operator==(const DiracMatrix& g) const {
  for (int i(0); i < 4; i++) {
    for (int j(0); j < 4; j++) {
      if (g.x[i][j] != x[i][j]) {
        return false;
      }
    }
  }
  return true;
}
/// Unequality check.
bool DiracMatrix::operator!=(const DiracMatrix& g) const {
  for (int i(0); i < 4; i++) {
    for (int j(0); j < 4; j++) {
      if (g.x[i][j] != x[i][j]) {
        return true;
      }
    }
  }
  return false;
}
/// Print the matrix to the given stream in the form '((x11, x12, x13), (...),
/// (...))'.
ostream& operator<<(ostream& out, const DiracMatrix& g) {
  out << "("
      << "(" << g.x[0][0] << ", " << g.x[0][1] << ", " << g.x[0][2] << ", "
      << g.x[0][3] << "), "
      << "(" << g.x[1][0] << ", " << g.x[1][1] << ", " << g.x[1][2] << ", "
      << g.x[1][3] << "), "
      << "(" << g.x[2][0] << ", " << g.x[2][1] << ", " << g.x[2][2] << ", "
      << g.x[2][3] << "), "
      << "(" << g.x[3][0] << ", " << g.x[3][1] << ", " << g.x[3][2] << ", "
      << g.x[3][3] << ")"
      << ")";
  return out;
}
/// Print to a stream in the form of a table.
void DiracMatrix::output(ostream& out, int fieldWidth) const {
  out << "(" << setw(fieldWidth) << x[0][0] << "   " << setw(fieldWidth)
      << x[0][1] << "   " << setw(fieldWidth) << x[0][2] << "   "
      << setw(fieldWidth) << x[0][3] << ")\n"
      << "(" << setw(fieldWidth) << x[1][0] << "   " << setw(fieldWidth)
      << x[1][1] << "   " << setw(fieldWidth) << x[1][2] << "   "
      << setw(fieldWidth) << x[1][3] << ")\n"
      << "(" << setw(fieldWidth) << x[2][0] << "   " << setw(fieldWidth)
      << x[2][1] << "   " << setw(fieldWidth) << x[2][2] << "   "
      << setw(fieldWidth) << x[2][3] << ")\n"
      << "(" << setw(fieldWidth) << x[3][0] << "   " << setw(fieldWidth)
      << x[3][1] << "   " << setw(fieldWidth) << x[3][2] << "   "
      << setw(fieldWidth) << x[3][3] << ")\n";
}
/// Trace.
dcomplex DiracMatrix::trace() const {
  dcomplex ret(0);
  for (int i(0); i < 4; i++) {
    ret += x[i][i];
  }
  return ret;
}

dcomplex trace(const DiracMatrix& g) {
  dcomplex ret(0);
  for (int i(0); i < 4; i++) {
    ret += g.x[i][i];
  }
  return ret;
}

DiracMatrix DiracMatrix::conj() const {
  DiracMatrix g;
  for (int i(0); i < 4; i++) {
    for (int j(0); j < 4; j++) {
      g.x[i][j] = std::conj(x[i][j]);
    }
  }
  return g;
}

DiracMatrix DiracMatrix::transpose() const {
  DiracMatrix g;
  for (int i(0); i < 4; i++) {
    for (int j(0); j < 4; j++) {
      g.x[i][j] = x[j][i];
    }
  }
  return g;
}

DiracMatrix DiracMatrix::trans_conj() const {
  DiracMatrix g;
  for (int i(0); i < 4; i++) {
    for (int j(0); j < 4; j++) {
      g.x[i][j] = std::conj(x[j][i]);
    }
  }
  return g;
}

const DiracMatrix gamma0_(tau_unit, tau_null, tau_null, -tau_unit);
const DiracMatrix gamma1_(tau_null, tau1, -tau1, tau_null);
const DiracMatrix gamma2_(tau_null, tau2, -tau2, tau_null);
const DiracMatrix gamma3_(tau_null, tau3, -tau3, tau_null);
const DiracMatrix gamma5_(tau_null, tau_unit, tau_unit, tau_null);
const DiracMatrix gamma_null(tau_null, tau_null, tau_null, tau_null);
const DiracMatrix gamma_unit(tau_unit, tau_null, tau_null, tau_unit);

const Gamma_mu& Gamma_mu::instance() { return gamma_; }

const DiracMatrix& Gamma_mu::operator()(int mu) const {
  assert(0 <= mu and mu < 4);
  if (mu == 0) {
    return gamma0_;
  }
  if (mu == 1) {
    return gamma1_;
  }
  if (mu == 2) {
    return gamma2_;
  }
  if (mu == 3) {
    return gamma3_;
  }
  return gamma_null;  // never get here
}

DiracMatrix Gamma_mu::operator()(const FourVector& p) const {
  DiracMatrix ret = gamma_null;
  for (int mu(0); mu < 4; mu++) {
    ret += sign_(mu) * gamma_(mu) * p(mu);
  }
  return ret;
}

Gamma_mu::Gamma_mu() {}

const Gamma_mu& gamma_ = Gamma_mu::instance();

const Gamma_mu Gamma_mu::gamma_;

const Sigma_mu_nu& Sigma_mu_nu::instance() { return sigma_; }

DiracMatrix Sigma_mu_nu::operator()(int mu, int nu) const {
  assert(0 <= mu and mu < 4);
  assert(0 <= nu and nu < 4);
  return i_ / 2. * (gamma_(mu) * gamma_(nu) - gamma_(nu) * gamma_(mu));
}
DiracMatrix Sigma_mu_nu::operator()(int mu, const FourVector& p) const {
  assert(0 <= mu and mu < 4);
  DiracMatrix ret = gamma_null;
  for (int nu(0); nu < 4; nu++) {
    ret += i_ / 2. * (gamma_(mu) * gamma_(nu) - gamma_(nu) * gamma_(mu)) *
           p(nu) * sign_(mu);
  }
  return ret;
}

DiracMatrix Sigma_mu_nu::operator()(const FourVector& p, int nu) const {
  assert(0 <= nu and nu < 4);
  DiracMatrix ret = gamma_null;
  for (int mu(0); mu < 4; mu++) {
    ret += i_ / 2. * (gamma_(mu) * gamma_(nu) - gamma_(nu) * gamma_(mu)) *
           p(mu) * sign_(mu);
  }
  return ret;
}

DiracMatrix Sigma_mu_nu::operator()(const FourVector& p1,
                                    const FourVector& p2) const {
  DiracMatrix ret = gamma_null;
  for (int mu(0); mu < 4; mu++) {
    for (int nu(0); nu < 4; nu++) {
      ret += i_ / 2. * (gamma_(mu) * gamma_(nu) - gamma_(nu) * gamma_(mu)) *
             p1(mu) * sign_(mu) * p2(nu) * sign_(nu);
    }
  }
  return ret;
}

Sigma_mu_nu::Sigma_mu_nu() {}

const Sigma_mu_nu& sigma_ = Sigma_mu_nu::instance();

const Sigma_mu_nu Sigma_mu_nu::sigma_;

}  // end of namespace Spinors
