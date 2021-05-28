#include <cmath>
#include <iomanip>
#include <cassert>
#include "Vectors.hpp"

using namespace std;

namespace Vectors {

  bool isIndex3(uint i) {
    return (1<=i and i<=3);
  }

  bool isIndex4(uint mu) {
    return (0<=mu and mu<=3);
  }

  ThreeVector::ThreeVector() {
    for (uint i(1); i<4; i++) {
      x[i] = 0.;
    }
  }

  ThreeVector::ThreeVector(const ThreeVector& other) {
    for (uint i(1); i<4; i++) {
      x[i] = other.x[i];
    }
  }
 
  ThreeVector::ThreeVector(double x1, double x2, double x3) {
    x[1] = x1;
    x[2] = x2;
    x[3] = x3;
  }

  ThreeVector::~ThreeVector() {}

  ThreeVector& ThreeVector::operator = (const ThreeVector& other) {
    if (&other != this) {
      for (uint i(1); i<4; i++) {
        x[i] = other.x[i];
      }
    }
    return *this;
  }

  double& ThreeVector::operator () (uint i) {
    assert(isIndex3(i));
    return x[i];
  }

  const double& ThreeVector::operator () (uint i) const {
    assert(isIndex3(i));
    return x[i];
  }

  const ThreeVector& ThreeVector::operator + () const {
    return *this;
  }

  const ThreeVector ThreeVector::operator - () const {
    return (-1)*(*this);
  }

  const ThreeVector ThreeVector::operator + (const ThreeVector& other) const {
    return ThreeVector(x[1]+other.x[1], x[2]+other.x[2], x[3]+other.x[3]);
  }

  const ThreeVector ThreeVector::operator - (const ThreeVector& other) const {
    return ThreeVector(x[1]-other.x[1], x[2]-other.x[2], x[3]-other.x[3]);
  }

  ThreeVector& ThreeVector::operator += (const ThreeVector& other) {
    for (uint i(1); i<4; i++) {
      x[i] += other.x[i];
    }
    return *this;
  }

  ThreeVector& ThreeVector::operator -= (const ThreeVector& other) {
    for (uint i(1); i<4; i++) {
      x[i] -= other.x[i];
    }
    return *this;
  }


  const double ThreeVector::operator * (const ThreeVector& other) const {
    return x[1]*other.x[1] + x[2]*other.x[2] + x[3]*other.x[3];
  }

  const ThreeVector ThreeVector::operator * (const double& c) const {
    return ThreeVector(x[1]*c, x[2]*c, x[3]*c);
  }

  const ThreeVector ThreeVector::operator / (const double& c) const {
    return ThreeVector(x[1]/c, x[2]/c, x[3]/c);
  }

  ThreeVector& ThreeVector::operator *= (const double& c) {
    for (uint i(1); i<4; i++) {
      x[i] *= c;
    }
    return *this;
  }

  ThreeVector& ThreeVector::operator /= (const double& c) {
    for (uint i(1); i<4; i++) {
      x[i] /= c;
    }
    return *this;
  }

  const ThreeVector operator * (const double& c, const ThreeVector& v) {
    return ThreeVector(c*v.x[1], c*v.x[2], c*v.x[3]);
  }

  ThreeVector crossprod(const ThreeVector& v1, const ThreeVector& v2) {
    return ThreeVector(v1.x[2]*v2.x[3] - v1.x[3]*v2.x[2],
                       v1.x[3]*v2.x[1] - v1.x[1]*v2.x[3],
                       v1.x[1]*v2.x[2] - v1.x[2]*v2.x[1]);
  }

  bool ThreeVector::operator == (const ThreeVector& other) const {
    return (x[1]==other.x[1]) && (x[2]==other.x[2]) && (x[3]==other.x[3]);
  }

  bool ThreeVector::operator != (const ThreeVector& other) const {
    return (x[1]!=other.x[1]) || (x[2]!=other.x[2]) || (x[3]!=other.x[3]);
  }

  ostream& operator << (ostream& out, const ThreeVector& v) {
    out << "(" << v.x[1] << ", " << v.x[2] << ", " << v.x[3] << ")";
    return out;
  }

  /**
     This is a very primitive implementation of operator>>, accepting only
     the form '[vx,vy,vz]' (no whitespaces etc.). [] are used because () is 
     interpreted by bash.Must be improved. (A separate parser class should 
     be used.)
  */
  istream& operator >> (istream& in, ThreeVector& v) {
    v = ThreeVector(0,0,0);
    ThreeVector x(0,0,0);
    char c;
    try {
      in.get(c);
      if (c!='[') { return in; }
      in >> x.x[1];
      in.get(c);
      if (c!=',') { return in; }
      in >> x.x[2];
      in.get(c);
      if (c!=',') { return in; }
      in >> x.x[3];
      in.get(c);
      if (c!=']') { return in; }
      v = x;
    } catch(...) {} 
    return in;
  }

  double ThreeVector::abs() const {
    return sqrt(x[1]*x[1] + x[2]*x[2] + x[3]*x[3]);
  }

  double ThreeVector::square() const {
    return x[1]*x[1] + x[2]*x[2] + x[3]*x[3];
  }

  bool ThreeVector::isNaN() const {
    return std::isnan(x[1]) or std::isnan(x[2]) or std::isnan(x[3]);
  }

  const ThreeVector ThreeVector::nullVector(0.,0.,0.);

  const ThreeVector ThreeVector::e[] = {
    ThreeVector(0,0,0),
    ThreeVector(1,0,0),
    ThreeVector(0,1,0),
    ThreeVector(0,0,1)
  };

  ThreeTensor::ThreeTensor() {
    for (uint i(1); i<4; i++) {
      for (uint j(1); j<4; j++) {
        x[i][j] = 0.;
      }
    }
  }

  ThreeTensor::ThreeTensor(const ThreeTensor& other) {
    for (uint i(1); i<4; i++) {
      for (uint j(1); j<4; j++) {
        x[i][j] = other.x[i][j];
      }
    }
  }
 
  ThreeTensor::ThreeTensor(double x11, double x12, double x13,
                           double x21, double x22, double x23,
                           double x31, double x32, double x33) {
    x[1][1] = x11;
    x[1][2] = x12;
    x[1][3] = x13;
    x[2][1] = x21;
    x[2][2] = x22;
    x[2][3] = x23;
    x[3][1] = x31;
    x[3][2] = x32;
    x[3][3] = x33;
  }

  ThreeTensor::ThreeTensor(double x1, double x2, double x3) {
    x[1][1] = x1;
    x[2][2] = x2;
    x[3][3] = x3;
    x[1][2] = 0.;
    x[1][3] = 0.;
    x[2][1] = 0.;
    x[2][3] = 0.;
    x[3][1] = 0.;
    x[3][2] = 0.;
  }

  ThreeTensor::~ThreeTensor() {}

  ThreeTensor& ThreeTensor::operator = (const ThreeTensor& other) {
    if (&other != this) {
      for (uint i(1); i<4; i++) {
        for (uint j(1); j<4; j++) {
          x[i][j] = other.x[i][j];
        }
      }
    }
    return *this;
  }

  double& ThreeTensor::operator () (uint i,uint j) {
    assert(isIndex3(i));
    assert(isIndex3(j));
    return x[i][j];
  }

  const double& ThreeTensor::operator () (uint i,uint j) const {
    assert(isIndex3(i));
    assert(isIndex3(j));
    return x[i][j];
  }

  const ThreeTensor& ThreeTensor::operator + () const {
    return *this;
  }

  const ThreeTensor ThreeTensor::operator - () const {
    return (-1)*(*this);
  }

  const ThreeTensor ThreeTensor::operator + (const ThreeTensor& other) const {
    return ThreeTensor(x[1][1]+other.x[1][1], x[1][2]+other.x[1][2], x[1][3]+other.x[1][3],
                       x[2][1]+other.x[2][1], x[2][2]+other.x[2][2], x[2][3]+other.x[2][3],
                       x[3][1]+other.x[3][1], x[3][2]+other.x[3][2], x[3][3]+other.x[3][3]);
  }

  const ThreeTensor ThreeTensor::operator - (const ThreeTensor& other) const {
    return ThreeTensor(x[1][1]-other.x[1][1], x[1][2]-other.x[1][2], x[1][3]-other.x[1][3],
                       x[2][1]-other.x[2][1], x[2][2]-other.x[2][2], x[2][3]-other.x[2][3],
                       x[3][1]-other.x[3][1], x[3][2]-other.x[3][2], x[3][3]-other.x[3][3]);
  }

  ThreeTensor& ThreeTensor::operator += (const ThreeTensor& other) {
    for (uint i(1); i<4; i++) {
      for (uint j(1); j<4; j++) {
        x[i][j] += other.x[i][j];
      }
    }
    return *this;
  }

  ThreeTensor& ThreeTensor::operator -= (const ThreeTensor& other) {
    for (uint i(1); i<4; i++) {
      for (uint j(1); j<4; j++) {
        x[i][j] -= other.x[i][j];
      }
    }
    return *this;
  }

  const ThreeTensor ThreeTensor::operator * (const ThreeTensor& other) const {
    return ThreeTensor(x[1][1]*other.x[1][1] + x[1][2]*other.x[2][1] + x[1][3]*other.x[3][1],
                       x[1][1]*other.x[1][2] + x[1][2]*other.x[2][2] + x[1][3]*other.x[3][2],
                       x[1][1]*other.x[1][3] + x[1][2]*other.x[2][3] + x[1][3]*other.x[3][3],
                       x[2][1]*other.x[1][1] + x[2][2]*other.x[2][1] + x[2][3]*other.x[3][1],
                       x[2][1]*other.x[1][2] + x[2][2]*other.x[2][2] + x[2][3]*other.x[3][2],
                       x[2][1]*other.x[1][3] + x[2][2]*other.x[2][3] + x[2][3]*other.x[3][3],
                       x[3][1]*other.x[1][1] + x[3][2]*other.x[2][1] + x[3][3]*other.x[3][1],
                       x[3][1]*other.x[1][2] + x[3][2]*other.x[2][2] + x[3][3]*other.x[3][2],
                       x[3][1]*other.x[1][3] + x[3][2]*other.x[2][3] + x[3][3]*other.x[3][3]);
  }

  const ThreeVector ThreeTensor::operator * (const ThreeVector& v) const {
    return ThreeVector(x[1][1]*v.x[1] + x[1][2]*v.x[2] + x[1][3]*v.x[3],
                       x[2][1]*v.x[1] + x[2][2]*v.x[2] + x[2][3]*v.x[3],
                       x[3][1]*v.x[1] + x[3][2]*v.x[2] + x[3][3]*v.x[3]);
  }

  const ThreeVector operator * (const ThreeVector& v, const ThreeTensor& t) {
    return ThreeVector(v.x[1]*t.x[1][1] + v.x[2]*t.x[2][1] + v.x[3]*t.x[3][1],
                       v.x[1]*t.x[1][2] + v.x[2]*t.x[2][2] + v.x[3]*t.x[3][2],
                       v.x[1]*t.x[1][3] + v.x[2]*t.x[2][3] + v.x[3]*t.x[3][3]);
  }

  const ThreeTensor diad(const ThreeVector& v, const ThreeVector& w) {
    return ThreeTensor(v.x[1]*w.x[1],v.x[1]*w.x[2],v.x[1]*w.x[3],
                       v.x[2]*w.x[1],v.x[2]*w.x[2],v.x[2]*w.x[3],
                       v.x[3]*w.x[1],v.x[3]*w.x[2],v.x[3]*w.x[3]);
  }

  const ThreeTensor ThreeTensor::operator * (const double& c) const {
    return ThreeTensor(x[1][1]*c, x[1][2]*c, x[1][3]*c,
                       x[2][1]*c, x[2][2]*c, x[2][3]*c,
                       x[3][1]*c, x[3][2]*c, x[3][3]*c);
  }

  const ThreeTensor ThreeTensor::operator / (const double& c) const {
    return ThreeTensor(x[1][1]/c, x[1][2]/c, x[1][3]/c,
                       x[2][1]/c, x[2][2]/c, x[2][3]/c,
                       x[3][1]/c, x[3][2]/c, x[3][3]/c);
  }

  ThreeTensor& ThreeTensor::operator *= (const double& c) {
    for (uint i(1); i<4; i++) {
      for (uint j(1); j<4; j++) {
        x[i][j] *= c;
      }
    }
    return *this;
  }

  ThreeTensor& ThreeTensor::operator /= (const double& c) {
    for (uint i(1); i<4; i++) {
      for (uint j(1); j<4; j++) {
        x[i][j] /= c;
      }
    }
    return *this;
  }

  const ThreeTensor operator * (const double& c, const ThreeTensor& t) {
    return ThreeTensor(c*t.x[1][1], c*t.x[1][2], c*t.x[1][3],
                       c*t.x[2][1], c*t.x[2][2], c*t.x[2][3],
                       c*t.x[3][1], c*t.x[3][2], c*t.x[3][3]);
  }

  bool ThreeTensor::operator == (const ThreeTensor& other) const {
    return ((x[1][1]==other.x[1][1]) && (x[1][2]==other.x[1][2]) && (x[1][3]==other.x[1][3]) &&
            (x[2][1]==other.x[2][1]) && (x[2][2]==other.x[2][2]) && (x[2][3]==other.x[2][3]) &&
            (x[3][1]==other.x[3][1]) && (x[3][2]==other.x[3][2]) && (x[3][3]==other.x[3][3]));
  }

  bool ThreeTensor::operator != (const ThreeTensor& other) const {
    return ((x[1][1]!=other.x[1][1]) || (x[1][2]!=other.x[1][2]) || (x[1][3]!=other.x[1][3]) ||
            (x[2][1]!=other.x[2][1]) || (x[2][2]!=other.x[2][2]) || (x[2][3]!=other.x[2][3]) ||
            (x[3][1]!=other.x[3][1]) || (x[3][2]!=other.x[3][2]) || (x[3][3]!=other.x[3][3]));
  }

  ostream& operator << (ostream& out, const ThreeTensor& t) {
      out << "((" << t.x[1][1] << ", " << t.x[1][2] << ", " << t.x[1][3] << "), "
          <<  "(" << t.x[2][1] << ", " << t.x[2][2] << ", " << t.x[2][3] << "), "
          <<  "(" << t.x[3][1] << ", " << t.x[3][2] << ", " << t.x[3][3] << "))";
    return out;
  }

  void ThreeTensor::output(ostream& out, uint fieldWidth) const {
    out << "\n/ " << setw(fieldWidth) << x[1][1] << "  "
        << setw(fieldWidth) << x[1][2] << "  "  << setw(fieldWidth) << x[1][3] << " \\\n"
        << "| " << setw(fieldWidth) << x[2][1] << "  "
        << setw(fieldWidth) << x[2][2] << "  "  << setw(fieldWidth) << x[2][3] << " |\n"
        << "\\ " << setw(fieldWidth) << x[3][1] << "  "
        << setw(fieldWidth) << x[3][2] << "  "  << setw(fieldWidth) << x[3][3] << " /" << endl;
  }

  double ThreeTensor::trace() const {
    return x[1][1] + x[2][2] + x[3][3];
  }

  ThreeTensor ThreeTensor::diagonalTensor(const ThreeVector& v) {
    return ThreeTensor(v.x[1],v.x[2],v.x[3]);
  }

  const ThreeTensor ThreeTensor::nullTensor(0.,0.,0.);

  const ThreeTensor ThreeTensor::unitTensor(1.,1.,1.);

  FourVector::FourVector() {
    for (uint i(0); i<4; i++) {
      x[i] = 0.;
    }
  }

  FourVector::FourVector(const FourVector& other) {
    for (uint i(0); i<4; i++) {
      x[i] = other.x[i];
    }
  }
 
  FourVector::FourVector(double x0, double x1, double x2, double x3) {
    x[0] = x0;
    x[1] = x1;
    x[2] = x2;
    x[3] = x3;
  }

  FourVector::FourVector(double x0, const ThreeVector& xx) {
    x[0] = x0;
    x[1] = xx(1);
    x[2] = xx(2);
    x[3] = xx(3);

  }

  FourVector::~FourVector() {}

  FourVector& FourVector::operator = (const FourVector& other) {
    if (&other != this) {
      for (uint i(0); i<4; i++) {
        x[i] = other.x[i];
      }
    }
    return *this;
  }

  double& FourVector::operator () (uint mu) {
    assert(isIndex4(mu));
    return x[mu];
  }

  const double& FourVector::operator () (uint mu) const {
    assert(isIndex4(mu));
    return x[mu];
  }

  const FourVector& FourVector::operator + () const {
    return *this;
  }

  const FourVector FourVector::operator - () const {
    return (-1)*(*this);
  }

  const FourVector FourVector::operator + (const FourVector& other) const {
    return FourVector(x[0]+other.x[0], x[1]+other.x[1], x[2]+other.x[2], x[3]+other.x[3]);
  }

  const FourVector FourVector::operator - (const FourVector& other) const {
    return FourVector(x[0]-other.x[0], x[1]-other.x[1], x[2]-other.x[2], x[3]-other.x[3]);
  }

  FourVector& FourVector::operator += (const FourVector& other) {
    for (uint i(0); i<4; i++) {
      x[i] += other.x[i];
    }
    return *this;
  }

  FourVector& FourVector::operator -= (const FourVector& other) {
    for (uint i(0); i<4; i++) {
      x[i] -= other.x[i];
    }
    return *this;
  }


  const double FourVector::operator * (const FourVector& other) const {
    return x[0]*other.x[0] - x[1]*other.x[1] - x[2]*other.x[2] - x[3]*other.x[3];
  }

  const FourVector FourVector::operator * (const double& c) const {
    return FourVector(x[0]*c, x[1]*c, x[2]*c, x[3]*c);
  }

  const FourVector FourVector::operator / (const double& c) const {
    return FourVector(x[0]/c, x[1]/c, x[2]/c, x[3]/c);
  }

  FourVector& FourVector::operator *= (const double& c) {
    for (uint i(0); i<4; i++) {
      x[i] *= c;
    }
    return *this;
  }

  FourVector& FourVector::operator /= (const double& c) {
    for (uint i(0); i<4; i++) {
      x[i] /= c;
    }
    return *this;
  }

  const FourVector operator * (const double& c, const FourVector& v) {
    return FourVector(c*v.x[0], c*v.x[1], c*v.x[2], c*v.x[3]);
  }

  bool FourVector::operator == (const FourVector& other) const {
    return (x[0]==other.x[0]) && (x[1]==other.x[1])
      && (x[2]==other.x[2]) && (x[3]==other.x[3]);
  }

  bool FourVector::operator != (const FourVector& other) const {
    return (x[0]!=other.x[0]) || (x[1]!=other.x[1])
      || (x[2]!=other.x[2]) || (x[3]!=other.x[3]);
  }

  ostream& operator << (ostream& out, const FourVector& v) {
    out << "(" << v.x[0] << ", " << v.x[1] << ", " << v.x[2] << ", " << v.x[3] << ")";
    return out;
  }

  istream& operator >> (istream& in, FourVector& v) {
    in >> v.x[0] >> v.x[1] >> v.x[2] >> v.x[3];
    return in;
  }

  double FourVector::square() const {
    return x[0]*x[0] - x[1]*x[1] - x[2]*x[2] - x[3]*x[3];
  }

  ThreeVector FourVector::spacial() const {
    return ThreeVector(x[1], x[2], x[3]);
  }

  const FourVector FourVector::reverse() const {
    return FourVector(x[0],-x[1],-x[2],-x[3]);
  }
  const FourVector FourVector::reflect() const {
    return FourVector(x[0],-x[1],-x[2],-x[3]);
  }

  FourVector reflect(const FourVector& v) {
    return FourVector(v.x[0],-v.x[1],-v.x[2],-v.x[3]);
  }

  bool FourVector::timelike() const {
    return (square() > 0);
  }

  bool FourVector::future() const {
    return ((square() > 0) and x[0]>0);
  }

  bool FourVector::past() const {
    return ((square() > 0) and x[0]<0);
  }

  bool FourVector::spacelike() const {
    return (square() < 0);
  }

  bool FourVector::isNaN() const {
    return std::isnan(x[0]) or std::isnan(x[1]) or std::isnan(x[2]) or std::isnan(x[3]);
  }

  const FourVector FourVector::nullVector(0.,0.,0.,0.);

  const FourVector FourVector::e[] = {
    FourVector(1,0,0,0),
    FourVector(0,1,0,0),
    FourVector(0,0,1,0),
    FourVector(0,0,0,1)
  };

  FourTensor::FourTensor() {
    for (uint i(0); i<4; i++) {
      for (uint j(0); j<4; j++) {
        x[i][j] = 0.;
      }
    }
  }

  FourTensor::FourTensor(const FourTensor& other) {
    for (uint i(0); i<4; i++) {
      for (uint j(0); j<4; j++) {
        x[i][j] = other.x[i][j];
      }
    }
  }
 
  FourTensor::FourTensor(double x00, double x01, double x02, double x03,
                         double x10, double x11, double x12, double x13,
                         double x20, double x21, double x22, double x23,
                         double x30, double x31, double x32, double x33) {
    x[0][0] = x00;
    x[0][1] = x01;
    x[0][2] = x02;
    x[0][3] = x03;
    x[1][0] = x10;
    x[1][1] = x11;
    x[1][2] = x12;
    x[1][3] = x13;
    x[2][0] = x20;
    x[2][1] = x21;
    x[2][2] = x22;
    x[2][3] = x23;
    x[3][0] = x30;
    x[3][1] = x31;
    x[3][2] = x32;
    x[3][3] = x33;
  }

    /// Static method to create FourTensor from mixed index components.
  FourTensor FourTensor::createUpDown(double x00, double x01, double x02, double x03,
                                      double x10, double x11, double x12, double x13,
                                      double x20, double x21, double x22, double x23,
                                      double x30, double x31, double x32, double x33) {
    return FourTensor(x00, -x01, -x02, -x03,
                      x10, -x11, -x12, -x13,
                      x20, -x21, -x22, -x23,
                      x30, -x31, -x32, -x33);
  }



  FourTensor::FourTensor(double x0, double x1, double x2, double x3) {
    x[0][0] = x0;
    x[0][1] = 0.;
    x[0][2] = 0.;
    x[0][3] = 0.;
    x[1][0] = 0.;
    x[1][1] = x1;
    x[1][2] = 0.;
    x[1][3] = 0.;
    x[2][0] = 0.;
    x[2][1] = 0.;
    x[2][2] = x2;
    x[2][3] = 0.;
    x[3][0] = 0.;
    x[3][1] = 0.;
    x[3][2] = 0.;
    x[3][3] = x3;
  }

  FourTensor::FourTensor(double x00, const ThreeVector& x0j,  const ThreeVector& xi0, const ThreeTensor& xij) {
    x[0][0] = x00;
    for (uint i(1); i<4; i++) {
      x[0][i] = x0j(i);
      x[i][0] = xi0(i);
      for (uint j(1); j<4; j++) {
        x[i][j] = xij(i,j);
      }
    }
  }

  FourTensor::~FourTensor() {}

  FourTensor& FourTensor::operator = (const FourTensor& other) {
    if (&other != this) {
      for (uint i(0); i<4; i++) {
        for (uint j(0); j<4; j++) {
          x[i][j] = other.x[i][j];
        }
      }
    }
    return *this;
  }

  double FourTensor::operator () (uint mu,uint nu) const {
    assert(isIndex4(mu));
    assert(isIndex4(nu));
    return x[mu][nu];
  }

  const FourVector FourTensor::operator () (uint mu) const {
    assert(isIndex4(mu));
    return FourVector(x[mu][0], x[mu][1], x[mu][2], x[mu][3]);
  }

  const FourTensor& FourTensor::operator + () const {
    return *this;
  }

  const FourTensor FourTensor::operator - () const {
    return (-1)*(*this);
  }

  const FourTensor FourTensor::operator + (const FourTensor& other) const {
    FourTensor ret;
    for (uint mu(0); mu<4; mu++) {
      for (uint nu(0); nu<4; nu++) {
        ret.x[mu][nu] = x[mu][nu] + other.x[mu][nu];
      }
    }
    return ret;
  }

  const FourTensor FourTensor::operator - (const FourTensor& other) const {
    FourTensor ret;
    for (uint mu(0); mu<4; mu++) {
      for (uint nu(0); nu<4; nu++) {
        ret.x[mu][nu] = x[mu][nu] - other.x[mu][nu];
      }
    }
    return ret;
  }

  FourTensor& FourTensor::operator += (const FourTensor& other) {
    for (uint i(0); i<4; i++) {
      for (uint j(0); j<4; j++) {
        x[i][j] += other.x[i][j];
      }
    }
    return *this;
  }

  FourTensor& FourTensor::operator -= (const FourTensor& other) {
    for (uint i(0); i<4; i++) {
      for (uint j(0); j<4; j++) {
        x[i][j] -= other.x[i][j];
      }
    }
    return *this;
  }

  const FourTensor FourTensor::operator * (const FourTensor& other) const {
    FourTensor ret;
    for (uint mu(0); mu<4; mu++) {
      for (uint nu(0); nu<4; nu++) {
        for (uint ro(0); ro<4; ro++) {
          ret.x[mu][nu] += sign_(ro) * x[mu][ro] * other.x[ro][nu];
        }
      }
    }
    return ret;
  }

  const FourVector FourTensor::operator * (const FourVector& v) const {
    FourVector ret;
    for (uint mu(0); mu<4; mu++) {
      for (uint nu(0); nu<4; nu++) {
        ret.x[mu] += x[mu][nu] * v.x[nu] * sign_(nu);
      }
    }
    return ret;
  }

  const FourVector operator * (const FourVector& v, const FourTensor& t) {
    FourVector ret;
    for (uint mu(0); mu<4; mu++) {
      for (uint nu(0); nu<4; nu++) {
        ret.x[mu] += sign_(nu) * v.x[nu] * t.x[nu][mu];
      }
    }
    return ret;
  }

  const FourTensor diad(const FourVector& v, const FourVector& w) {
    FourTensor ret;
    for (uint mu(0); mu<4; mu++) {
      for (uint nu(0); nu<4; nu++) {
        ret.x[mu][nu] = v.x[mu] * w.x[nu];
      }
    }
    return ret;
  }

  const FourTensor FourTensor::operator * (const double& c) const {
    FourTensor ret;
    for (uint mu(0); mu<4; mu++) {
      for (uint nu(0); nu<4; nu++) {
        ret.x[mu][nu] = x[mu][nu] * c;
      }
    }
    return ret;
  }

  const FourTensor FourTensor::operator / (const double& c) const {
    FourTensor ret;
    for (uint mu(0); mu<4; mu++) {
      for (uint nu(0); nu<4; nu++) {
        ret.x[mu][nu] = x[mu][nu] / c;
      }
    }
    return ret;
  }

  FourTensor& FourTensor::operator *= (const double& c) {
    for (uint i(0); i<4; i++) {
      for (uint j(0); j<4; j++) {
        x[i][j] *= c;
      }
    }
    return *this;
  }

  FourTensor& FourTensor::operator /= (const double& c) {
    for (uint i(0); i<4; i++) {
      for (uint j(0); j<4; j++) {
        x[i][j] /= c;
      }
    }
    return *this;
  }

  const FourTensor operator * (const double& c, const FourTensor& t) {
    FourTensor ret;
    for (uint mu(0); mu<4; mu++) {
      for (uint nu(0); nu<4; nu++) {
        ret.x[mu][nu] = t.x[mu][nu] * c;
      }
    }
    return ret;
  }

  bool FourTensor::operator == (const FourTensor& other) const {
    for (uint mu(0); mu<4; mu++) {
      for (uint nu(0); nu<4; nu++) {
        if (x[mu][nu] != other.x[mu][nu]) { return false; }
      }
    }
    return true;
  }

  bool FourTensor::operator != (const FourTensor& other) const {
    return not(*this == other);
  }

  ostream& operator << (ostream& out, const FourTensor& t) {
    out << "((" << t.x[0][0] << ", " << t.x[0][1] << ", " << t.x[0][2] << ", " << t.x[0][3] << "), "
        << " (" << t.x[1][0] << ", " << t.x[1][1] << ", " << t.x[1][2] << ", " << t.x[1][3] << "), "
        <<  "(" << t.x[2][0] << ", " << t.x[2][1] << ", " << t.x[2][2] << ", " << t.x[2][3] << "), "
        <<  "(" << t.x[3][0] << ", " << t.x[3][1] << ", " << t.x[3][2] << ", " << t.x[3][3] << "))";
    return out;
  }

  void FourTensor::output(ostream& out, uint fieldWidth) const {
    out << "\n/ " << setw(fieldWidth) << x[0][0] << "  "  << setw(fieldWidth) << x[0][1] << "  "
        << setw(fieldWidth) << x[0][2] << "  "  << setw(fieldWidth) << x[0][3] << " \\\n"
        << "| " << setw(fieldWidth) << x[1][0] << "  "  << setw(fieldWidth) << x[1][1] << "  "
        << setw(fieldWidth) << x[1][2] << "  "  << setw(fieldWidth) << x[1][3] << " |\n"
        << "| " << setw(fieldWidth) << x[2][0] << "  "  << setw(fieldWidth) << x[2][1] << "  "
        << setw(fieldWidth) << x[2][2] << "  "  << setw(fieldWidth) << x[2][3] << " |\n"
        << "\\ " << setw(fieldWidth) << x[3][0] << "  "  << setw(fieldWidth) << x[3][1] << "  "
        << setw(fieldWidth) << x[3][2] << "  "  << setw(fieldWidth) << x[3][3] << " /" << endl;
  }

  double FourTensor::trace() const {
    return x[0][0] + x[1][1] + x[2][2] + x[3][3];
  }

  const FourTensor FourTensor::nullTensor(0.,0.,0.,0.);

  const FourTensor FourTensor::unitTensor(1.,1.,1.,1.);

  LorentzTransformation::LorentzTransformation(const FourVector& a) : a(a) {
    const double msqr(a*a);
    if (msqr>0) { // a is timelike
      if (a.spacial() == ThreeVector::nullVector) {
        L = FourTensor::unitTensor;
      } else {
        const double m = sqrt(msqr);
        const ThreeVector aa = a.spacial();
        const ThreeTensor Lspacial = ThreeTensor::unitTensor + (a(0)/m-1.)*diad(aa,aa)/(aa*aa);
        L = FourTensor(a(0)/m,aa/m,aa/m,Lspacial);
        /*
        const FourVector u = a/m;
        const ThreeVector uu = u.spacial();
        const double gamma = u[0];
        const ThreeTensor Lspacial = ThreeTensor::unitTensor + diad(uu,uu)/(gamma+1.);
        L = FourTensor(gamma,uu,-uu,-Lspacial);
        */
      }
    } else if (msqr<0) { // a is spacelike
      if (a(0) == 0) {
        L = FourTensor::unitTensor;
      } else {
        const double m = sqrt(-msqr);
        const ThreeVector aa = a.spacial();
        const double aaabs = aa.abs();
        const ThreeTensor Lspacial = ThreeTensor::unitTensor + (aaabs/m-1.)*diad(aa,aa)/(aa*aa);
        L = FourTensor(aaabs/m,a(0)/(m*aaabs)*aa,a(0)/(m*aaabs)*aa,Lspacial);
      }
    } else {
      /*
        cerr << "LorentzTransformation::constr() - a is light-like or NaN" << endl;
        PR(a);
      //*/
      L = FourTensor::unitTensor;
    }
  }

  FourVector LorentzTransformation::operator()(const FourVector& x) const {
    return L*x;
  }

  FourTensor LorentzTransformation::operator()(const FourTensor& T) const {
    return L*T*L;
  }


  LorentzTransformation LorentzTransformation::inverse() const {
    return LorentzTransformation(FourVector(a.reverse()));
  }

  ostream& operator<<(ostream& out, const LorentzTransformation& lt) {
    out << lt.L;
    return out;
  }

  LorentzTransformation LorentzTransformation::toRestFrame(const FourVector& a) {
    return LorentzTransformation(a);
  }

  LorentzTransformation LorentzTransformation::fromRestFrame(const FourVector& a) {
    return LorentzTransformation(a.reverse());
  }

  double sign_(uint mu) {
    assert(mu>=0 and mu<4);
    if (mu==0) { return 1; }
    return -1;
  }

  double sign_(const vector<uint>& mu) {
    double ret(1);
    for (uint i(0); i<mu.size(); i++) {
      assert(0<=mu[i] and mu[i]<4);
      if (mu[i] > 0) { ret *= -1; }
    }
    return ret;
  }

  const FourTensor g_(1,-1,-1,-1);

  int EpsilonTensor(uint mu, uint nu, uint ro, uint si) {
    uint al[4] = {mu,nu,ro,si};
    for (uint i(0); i<3; i++) {
      for (uint j(i+1); j<4; j++) {
        if (al[i]==al[j]) { return 0; }
      }
    }
    int ret = 1;
    for (uint max(3); max>0; max--) {
      for (uint j(0); j<max; j++) {
        if (al[j]>al[j+1]) {
          uint h = al[j];
          al[j] = al[j+1];
          al[j+1] = h;
          ret *= -1;
        }
      }
    }
    return ret;
  }

} // end of namespace Vectors
