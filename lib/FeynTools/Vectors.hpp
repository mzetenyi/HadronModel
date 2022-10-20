#ifndef VECTORS_HPP
#define VECTORS_HPP

#include <cmath>
#include <iostream>
#include <vector>

namespace Vectors {

bool isIndex3(uint);
bool isIndex4(uint);

/**
   Forward declarations.
*/
class ThreeVector;
class ThreeTensor;
class FourVector;
class FourTensor;

/**
   Represent a threevector.
   Has the usual functionality (addition, multiplication with scalar,
   inner product, etc.)
*/
class ThreeVector {
 public:
  friend class ThreeTensor;
  /// Create a nullvector.
  ThreeVector();
  /// Copy constructor.
  ThreeVector(const ThreeVector&);
  /// Create a vector with the given components.
  ThreeVector(double x1, double x2, double x3);
  /// Destructor. Does nothing.
  ~ThreeVector();
  /// Assignment operator.
  ThreeVector& operator=(const ThreeVector&);
  /// Return vector component.
  double& operator()(uint);
  const double& operator()(uint) const;
  /// Positive sign.
  const ThreeVector& operator+() const;
  /// Negative sign.
  ThreeVector operator-() const;
  /// Vector addition.
  ThreeVector operator+(const ThreeVector&) const;
  /// Vector subtraction.
  ThreeVector operator-(const ThreeVector&) const;
  /// Increment operator.
  ThreeVector& operator+=(const ThreeVector&);
  /// Decrement operator.
  ThreeVector& operator-=(const ThreeVector&);
  /// Inner product.
  double operator*(const ThreeVector&) const;
  /// Multiplication with scalar (double).
  ThreeVector operator*(const double&) const;
  /// Division by scalar.
  ThreeVector operator/(const double&) const;
  ThreeVector& operator*=(const double&);
  ThreeVector& operator/=(const double&);
  /// Friend operator for scalar*vector.
  friend ThreeVector operator*(const double&, const ThreeVector&);
  /// Vector*tensor product.
  friend ThreeVector operator*(const ThreeVector&, const ThreeTensor&);
  /// Vectorial product
  friend ThreeVector crossprod(const ThreeVector&, const ThreeVector&);
  /// Diadic product.
  friend ThreeTensor diad(const ThreeVector&, const ThreeVector&);
  /// Equality check.
  bool operator==(const ThreeVector&) const;
  /// Unequality check.
  bool operator!=(const ThreeVector&) const;
  /// Print the vector to the given stream in the form '(x,y)'.
  friend std::ostream& operator<<(std::ostream&, const ThreeVector&);
  /// Read the vector from the given stream. It must be in the form '(x,y)'.
  friend std::istream& operator>>(std::istream&, ThreeVector&);
  /// Absolute value.
  double abs() const;
  /// Square of the vector.
  double square() const;
  /// Decide if any component is NaN.
  bool isNaN() const;
  /// Null vector.
  static const ThreeVector nullVector;
  /// Basis vectors.
  static const ThreeVector e[];

 private:
  /// Variables to store the vector components. (x[0] is not used)
  double x[4];
};

/**
   Represent a threetensor.
   Has the usual functionality (addition, multiplication with scalar,
   tensor product, etc.)
*/
class ThreeTensor {
 public:
  friend class ThreeVector;
  /// Create a nulltensor.
  ThreeTensor();
  /// Copy constructor.
  ThreeTensor(const ThreeTensor&);
  /// Create a tensor from components.
  ThreeTensor(double, double, double, double, double, double, double, double,
              double);
  /// Create a diagonal tensor with the given components.
  ThreeTensor(double, double, double);
  /// Destructor. Does nothing.
  ~ThreeTensor();
  /// Assignment operator.
  ThreeTensor& operator=(const ThreeTensor&);
  /// Return tensor component.
  double& operator()(uint, uint);
  const double& operator()(uint, uint) const;
  /// Positive sign.
  const ThreeTensor& operator+() const;
  /// Negative sign.
  ThreeTensor operator-() const;
  /// Tensor addition.
  ThreeTensor operator+(const ThreeTensor&) const;
  /// Tensor subtraction.
  ThreeTensor operator-(const ThreeTensor&) const;
  /// Increment operator.
  ThreeTensor& operator+=(const ThreeTensor&);
  /// Decrement operator.
  ThreeTensor& operator-=(const ThreeTensor&);
  /// Tensor product.
  ThreeTensor operator*(const ThreeTensor&) const;
  /// Tensor*vector product.
  ThreeVector operator*(const ThreeVector&) const;
  /// Vector*tensor product.
  friend ThreeVector operator*(const ThreeVector&, const ThreeTensor&);
  /// Multiplication with scalar (double).
  ThreeTensor operator*(const double&) const;
  /// Division by scalar.
  ThreeTensor operator/(const double&) const;
  ThreeTensor& operator*=(const double&);
  ThreeTensor& operator/=(const double&);
  /// Friend operator for scalar*tensor.
  friend ThreeTensor operator*(const double&, const ThreeTensor&);
  /// Equality check.
  bool operator==(const ThreeTensor&) const;
  /// Unequality check.
  bool operator!=(const ThreeTensor&) const;
  /// Print the tensor to the given stream in the form
  /// '((x11,x12,x13),(...),(...))'.
  friend std::ostream& operator<<(std::ostream&, const ThreeTensor&);
  /// Print to a stream in the form of a table.
  void output(std::ostream&, uint fieldWidth = 10) const;
  /// Trace.
  double trace() const;
  /// create a diagonal tensor with the components of the given vector
  static ThreeTensor diagonalTensor(const ThreeVector&);
  /// Null tensor.
  static const ThreeTensor nullTensor;
  /// Unit tensor.
  static const ThreeTensor unitTensor;

 private:
  /// Variables to store the tensor components. (x[0] is not used)
  double x[4][4];
};

/**
   Represent a fourvector.
   Has the usual functionality (addition, multiplication with scalar,
   inner product, etc.)

   \todo spatial() could return a reference (FourVector and
   ThreeVector have exactly equal representations.
*/
class FourVector {
 public:
  /// Create a nullvector.
  FourVector();
  /// Copy constructor.
  FourVector(const FourVector&);
  /// Create a vector with the given components.
  FourVector(double x0, double x1, double x2, double x3);
  /// Create a vector from a three scalar and a threevector.
  FourVector(double x0, const ThreeVector& xx);
  /// Destructor. Does nothing.
  ~FourVector();
  /// Assignment operator.
  FourVector& operator=(const FourVector&);
  /// Return contravariant vector component.
  double& operator()(uint);
  /// Return contravariant vector component.
  const double& operator()(uint) const;
  /// Positive sign.
  const FourVector& operator+() const;
  /// Negative sign.
  FourVector operator-() const;
  /// Vector addition.
  FourVector operator+(const FourVector&) const;
  /// Vector subtraction.
  FourVector operator-(const FourVector&) const;
  /// Increment operator.
  FourVector& operator+=(const FourVector&);
  /// Decrement operator.
  FourVector& operator-=(const FourVector&);
  /// Inner product.
  double operator*(const FourVector&) const;
  /// Vector*tensor product.
  friend FourVector operator*(const FourVector&, const FourTensor&);
  /// Diadic product.
  friend FourTensor diad(const FourVector&, const FourVector&);
  /// Multiplication with scalar (double).
  FourVector operator*(const double&) const;
  /// Division by scalar.
  FourVector operator/(const double&) const;
  FourVector& operator*=(const double&);
  FourVector& operator/=(const double&);
  /// Friend operator for scalar*vector.
  friend FourVector operator*(const double&, const FourVector&);
  /// Equality check.
  bool operator==(const FourVector&) const;
  /// Unequality check.
  bool operator!=(const FourVector&) const;
  /// Print the vector to the given stream.
  friend std::ostream& operator<<(std::ostream&, const FourVector&);
  /// Read in components of the vector from the given stream.
  friend std::istream& operator>>(std::istream&, FourVector&);
  /// Square of the vector.
  double square() const;
  /// Return the spacelike part as a ThreeVector.
  ThreeVector spacial() const;
  /// Reverse the spacial part of the FourVector. (Same as reflect())
  FourVector reverse() const;
  /// Perform a space reflection on the FourVector.
  FourVector reflect() const;
  /// Perform a space reflection on the FourVector.
  friend FourVector reflect(const FourVector&);
  /// Decide if the FourVector is timelike.
  bool timelike() const;
  /// Decide if the FourVector is timelike and future pointing.
  bool future() const;
  /// Decide if the FourVector is timelike and past pointing.
  bool past() const;
  /// Decide if the FourVector is spacelike.
  bool spacelike() const;
  /// Decide if any component is NaN.
  bool isNaN() const;
  /// Null vector.
  static const FourVector nullVector;
  /// Basis vectors.
  static const FourVector e[];
  /// Declare FourTensor as friend.
  friend class FourTensor;

 private:
  /// Variables to store the vector components.
  double x[4];
};

/**
   Represent a fourtensor.
   Has the usual functionality (addition, multiplication with scalar,
   inner product, etc.)
   The components are stored in the array x[0..4][0..4].
   All indices are upper.
*/
class FourTensor {
 public:
  /// Create a nulltensor.
  FourTensor();
  /// Copy constructor.
  FourTensor(const FourTensor&);
  /// Create a tensor from components. (Indices are upper.)
  FourTensor(double, double, double, double, double, double, double, double,
             double, double, double, double, double, double, double, double);
  /// Static method to create FourTensor from mixed index components.
  static FourTensor createUpDown(double, double, double, double, double, double,
                                 double, double, double, double, double, double,
                                 double, double, double, double);
  /// Create a diagonal tensor with the given components. (Indices are upper.)
  FourTensor(double x0, double x1, double x2, double x3);
  /// Create a Fourtensor from blocks. (Indices are upper.)
  FourTensor(double x00, const ThreeVector& x0j, const ThreeVector& xi0,
             const ThreeTensor& xij);
  /// Destructor. Does nothing.
  ~FourTensor();
  /// Assignment operator.
  FourTensor& operator=(const FourTensor&);
  /// Return upper index tensor component.
  double operator()(uint, uint) const;
  /// Return upper index components in a FourVector.
  FourVector operator()(uint) const;
  /// Positive sign.
  const FourTensor& operator+() const;
  /// Negative sign.
  FourTensor operator-() const;
  /// Tensor addition.
  FourTensor operator+(const FourTensor&) const;
  /// Tensor subtraction.
  FourTensor operator-(const FourTensor&) const;
  /// Increment operator.
  FourTensor& operator+=(const FourTensor&);
  /// Decrement operator.
  FourTensor& operator-=(const FourTensor&);
  /// Tensor product.
  FourTensor operator*(const FourTensor&) const;
  /// Tensor*vector product.
  FourVector operator*(const FourVector&) const;
  /// Vector*tensor product.
  friend FourVector operator*(const FourVector&, const FourTensor&);
  /// Diadic product.
  friend FourTensor diad(const FourVector&, const FourVector&);
  /// Multiplication with scalar (double).
  FourTensor operator*(const double&) const;
  /// Division by scalar.
  FourTensor operator/(const double&) const;
  FourTensor& operator*=(const double&);
  FourTensor& operator/=(const double&);
  /// Friend operator for scalar*tensor.
  friend FourTensor operator*(const double&, const FourTensor&);
  /// Equality check.
  bool operator==(const FourTensor&) const;
  /// Unequality check.
  bool operator!=(const FourTensor&) const;
  /// Print the tensor to the given stream in the form
  /// '(x00,x01,x02,x03),(..),(..),(..))'.
  friend std::ostream& operator<<(std::ostream&, const FourTensor&);
  /// Print to a stream in the form of a table.
  void output(std::ostream&, uint fieldWidth = 10) const;
  /// Trace.
  double trace() const;
  /// Null tensor.
  static const FourTensor nullTensor;
  /// Unit tensor.
  static const FourTensor unitTensor;
  /// Declare FourVector as friend.
  friend class FourVector;

 private:
  /// Variables to store the tensor components.
  double x[4][4];
};

/**
   Represents a Lorentz transformation that transforms fourvectors into
   the rest frame of the fourvector a.
*/
class LorentzTransformation {
 public:
  /**
     Constructor that creates a LorentzTransformation transforming to
     the rest frame of a.
  */
  LorentzTransformation(const FourVector& a);
  /// Transform the fourvector x to the rest frame of a.
  FourVector operator()(const FourVector& x) const;
  /// Transform the fourtensor T to the rest frame of a.
  FourTensor operator()(const FourTensor& T) const;
  /**
     Transform the quantity x to the rest frame of a. This operator
     allows the use of operator() to perform the Lorentz
     transfromation on an object of any type T, that has a member
     function transform(const LorentzTransformation&) - this member
     implements the transformation.
  */
  template <typename T>
  T operator()(const T& x) const {
    return x.transform(*this);
  }
  /// The inverse Lorentz transformation.
  LorentzTransformation inverse() const;

  /**
     Print components of the Lorentz tranformation matrix.
  */
  friend std::ostream& operator<<(std::ostream&, const LorentzTransformation&);
  /**
     Generate a LorentzTransformation that transforms to the rest
     frame of the given FourVector a.
  */
  static LorentzTransformation toRestFrame(const FourVector& a);
  /**
     Generate a LorentzTransformation that transforms a FourVector at
     rest to a FourVector parallel to the given FourVector a.
  */
  static LorentzTransformation fromRestFrame(const FourVector& a);

 private:
  const FourVector a;
  FourTensor L;
};

/**
   Signature of the metric.
*/
double sign_(uint mu);

/**
   Product of signs corresponding to the indices given.
*/
double sign_(const std::vector<uint>& mu);

/**
   Metric tensor.
*/
extern const FourTensor g_;

/**
   Totally antisymmetric epsilon tensor. EpsilonTensor(0,1,2,3) = 1,
   i.e. indices are upper.
*/
int EpsilonTensor(uint mu, uint nu, uint ro, uint si);

}  // namespace Vectors

#endif
