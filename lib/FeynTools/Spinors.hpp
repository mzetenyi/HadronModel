#ifndef SPINORS_HPP
#define SPINORS_HPP

#include <iostream>
#include <cmath>
#include "utils.hpp"
#include "Vectors.hpp"
#include "halfint.hpp"

namespace Spinors {

  class PauliSpinor;
  class AdPauliSpinor;
  class PauliMatrix;
  class DiracSpinor;
  class AdDiracSpinor;
  class DiracMatrix;

  /**
     Represent a two-component (Pauli) spinor.
     Has the usual functionality (addition, multiplication with scalar, 
     inner product, etc.) 
  */
  class PauliSpinor {
  public:
    /// Create a nullspinor.
    PauliSpinor();
    /// Copy constructor.
    PauliSpinor(const PauliSpinor&);
    /// Create a spinor with the given components.
    PauliSpinor(dcomplex,dcomplex);
    /// Destructor. Does nothing.
    ~PauliSpinor();
    /// Assignment operator.
    PauliSpinor& operator = (const PauliSpinor&);
    /// Adjoint
    const AdPauliSpinor adj() const;
    friend const AdPauliSpinor adj(const PauliSpinor&);
    /// Positive sign.
    const PauliSpinor& operator + () const;
    /// Negative sign.
    const PauliSpinor operator - () const;
    /// Spinor addition.
    const PauliSpinor operator + (const PauliSpinor&) const;
    /// Spinor subtraction.
    const PauliSpinor operator - (const PauliSpinor&) const;
    /// Increment operator.
    PauliSpinor& operator += (const PauliSpinor&);
    /// Decrement operator.
    PauliSpinor& operator -= (const PauliSpinor&);
    /// Multiplication with scalar (dcomplex).
    const PauliSpinor operator * (const dcomplex&) const;
    /// Division by scalar.
    const PauliSpinor operator / (const dcomplex&) const;
    PauliSpinor& operator *= (const dcomplex&);
    PauliSpinor& operator /= (const dcomplex&);
    /// Friend operator for scalar*spinor.
    friend const PauliSpinor operator * (const dcomplex&, const PauliSpinor&);
    /// Matrix*spinor product.
    friend const PauliSpinor operator * (const PauliMatrix&, const PauliSpinor&);
    /// Scalar product.
    friend const dcomplex operator * (const AdPauliSpinor&, const PauliSpinor&);
    /// Diadic product.
    friend const PauliMatrix operator * (const PauliSpinor&, const AdPauliSpinor&);
    /// Equality check.
    bool operator == (const PauliSpinor&) const;
    /// Unequality check.
    bool operator != (const PauliSpinor&) const;
    /// Print the spinor to the given stream in the form '(x,y)'.
    friend std::ostream& operator << (std::ostream&, const PauliSpinor&);
    /// Read the spinor from the given stream. It must be in the form '(x,y)'.
    friend std::istream& operator >> (std::istream&, PauliSpinor&);
    /// Decide if any component is NaN.
    bool isnan() const;
    /// friend classes for easier implementation.
    friend class DiracSpinor;
    friend class AdDiracSpinor;
    friend class DiracMatrix;
  private:
    /// Variables to store the spinor components.
    dcomplex x1,x2;
  };

  extern const PauliSpinor Khi_plus;
  extern const PauliSpinor Khi_minus;
  extern const PauliSpinor Khi_null;

  const PauliSpinor& Khi(halfint);

  /**
     Represent a two-component adjoint (Pauli) spinor.
     Has the usual functionality (addition, multiplication with scalar, 
     inner product, etc.) 
  */
  class AdPauliSpinor {
  public:
    /// Create a nullspinor.
    AdPauliSpinor();
    /// Copy constructor.
    AdPauliSpinor(const AdPauliSpinor&);
    /// Create a spinor with the given components.
    AdPauliSpinor(dcomplex,dcomplex);
    /// Destructor. Does nothing.
    ~AdPauliSpinor();
    /// Assignment operator.
    AdPauliSpinor& operator = (const AdPauliSpinor&);
    /// Adjoint
    const PauliSpinor adj() const;
    friend const PauliSpinor adj(const AdPauliSpinor&);
    /// Positive sign.
    const AdPauliSpinor& operator + () const;
    /// Negative sign.
    const AdPauliSpinor operator - () const;
    /// Spinor addition.
    const AdPauliSpinor operator + (const AdPauliSpinor&) const;
    /// Spinor subtraction.
    const AdPauliSpinor operator - (const AdPauliSpinor&) const;
    /// Increment operator.
    AdPauliSpinor& operator += (const AdPauliSpinor&);
    /// Decrement operator.
    AdPauliSpinor& operator -= (const AdPauliSpinor&);
    /// Multiplication with scalar (dcomplex).
    const AdPauliSpinor operator * (const dcomplex&) const;
    /// Division by scalar.
    const AdPauliSpinor operator / (const dcomplex&) const;
    AdPauliSpinor& operator *= (const dcomplex&);
    AdPauliSpinor& operator /= (const dcomplex&);
    /// Friend operator for scalar*spinor.
    friend const AdPauliSpinor operator * (const dcomplex&, const AdPauliSpinor&);
    /// Spinor*matrix product.
    friend const AdPauliSpinor operator * (const AdPauliSpinor&, const PauliMatrix&);
    /// Scalar product.
    friend const dcomplex operator * (const AdPauliSpinor&, const PauliSpinor&);
    /// Diadic product.
    friend const PauliMatrix operator * (const PauliSpinor&, const AdPauliSpinor&);
    /// Equality check.
    bool operator == (const AdPauliSpinor&) const;
    /// Unequality check.
    bool operator != (const AdPauliSpinor&) const;
    /// Print the spinor to the given stream in the form '(x,y)'.
    friend std::ostream& operator << (std::ostream&, const AdPauliSpinor&);
    /// Read the spinor from the given stream. It must be in the form '(x,y)'.
    friend std::istream& operator >> (std::istream&, AdPauliSpinor&);
    /// Decide if any component is NaN.
    bool isnan() const;
    /// friend classes for easier implementation.
    friend class DiracSpinor;
    friend class AdDiracSpinor;
    friend class DiracMatrix;
  private:
    /// Variables to store the spinor components.
    dcomplex x1,x2;
  };


  /**
     Represent a Pauli matrix.
  */
  class PauliMatrix {
  public:
    /// Create a nullmatrix.
    PauliMatrix();
    /// Copy constructor.
    PauliMatrix(const PauliMatrix&);
    /// Create a Pauli matrix from components.
    PauliMatrix(dcomplex,dcomplex,dcomplex,dcomplex);
    /// Destructor. Does nothing.
    ~PauliMatrix();
    /// Assignment operator.
    PauliMatrix& operator = (const PauliMatrix&);
    /// Adjoint
    const PauliMatrix adj() const;
    friend const PauliMatrix adj(const PauliMatrix&);
    /// Positive sign.
    const PauliMatrix& operator + () const;
    /// Negative sign.
    const PauliMatrix operator - () const;
    /// Matrix addition.
    const PauliMatrix operator + (const PauliMatrix&) const;
    /// Matrix subtraction.
    const PauliMatrix operator - (const PauliMatrix&) const;
    /// Increment operator.
    PauliMatrix& operator += (const PauliMatrix&);
    /// Decrement operator.
    PauliMatrix& operator -= (const PauliMatrix&);
    /// Matrix product.
    friend const PauliMatrix operator * (const PauliMatrix&, const PauliMatrix&);
    /// Matrix*spinor product.
    friend const PauliSpinor operator * (const PauliMatrix&, const PauliSpinor&);
    /// Spinor*matrix product.
    friend const AdPauliSpinor operator * (const AdPauliSpinor&, const PauliMatrix&);
    /// Multiplication with scalar (dcomplex).
    const PauliMatrix operator * (const dcomplex&) const;
    /// Division by scalar.
    const PauliMatrix operator / (dcomplex) const;
    PauliMatrix& operator *= (const dcomplex&);
    PauliMatrix& operator /= (const dcomplex&);
    /// Friend operator for scalar*matrix.
    friend const PauliMatrix operator * (const dcomplex&, const PauliMatrix&);
    /// Equality check.
    bool operator == (const PauliMatrix&) const;
    /// Unequality check.
    bool operator != (const PauliMatrix&) const;
    /// Print the matrix to the given stream in the form '((x11,x12,x13),(...),(...))'.
    friend std::ostream& operator << (std::ostream&, const PauliMatrix&);
    /// Print to a stream in the form of a table.
    void output(std::ostream&, int fieldWidth=10) const;
    /// Trace.
    dcomplex trace() const;
    /// friend classes for easier implementation.
    friend class DiracSpinor;
    friend class AdDiracSpinor;
    friend class DiracMatrix;
  private:
    /// Variables to store the matrix components.
    dcomplex x11,x12,x21,x22;
    const PauliMatrix conj() const;
    const PauliMatrix transpose() const;
    const PauliMatrix tran_conj() const;
  };

  extern const PauliMatrix tau1;
  extern const PauliMatrix tau2;
  extern const PauliMatrix tau3;
  extern const PauliMatrix tau_null;
  extern const PauliMatrix tau_unit;

  class Tau_i {
  public:
    static const Tau_i& instance();
    const PauliMatrix operator()(int) const;
    const PauliMatrix& operator[](int) const;
  private:
    static const Tau_i tau_;
    Tau_i();
    Tau_i(const Tau_i&);
  };

  extern const Tau_i& tau_;
  /**
     Represent a four-component (Dirac) spinor.
     Has the usual functionality (addition, multiplication with scalar, 
     inner product, adjoint, etc.) 
  */
  class DiracSpinor {
  public:
    /// Create a nullspinor.
    DiracSpinor();
    /// Copy constructor.
    DiracSpinor(const DiracSpinor&);
    /// Create a spinor with the given components.
    DiracSpinor(dcomplex,dcomplex,dcomplex,dcomplex);
    /// Create a spinor from two Pauli spinors.
    DiracSpinor(const PauliSpinor&, const PauliSpinor&);
    /// Destructor. Does nothing.
    ~DiracSpinor();
    /// Assignment operator.
    DiracSpinor& operator = (const DiracSpinor&);
    /// Positive sign.
    const DiracSpinor& operator + () const;
    /// Negative sign.
    const DiracSpinor operator - () const;
    /// Spinor addition.
    const DiracSpinor operator + (const DiracSpinor&) const;
    /// Spinor subtraction.
    const DiracSpinor operator - (const DiracSpinor&) const;
    /// Increment operator.
    DiracSpinor& operator += (const DiracSpinor&);
    /// Decrement operator.
    DiracSpinor& operator -= (const DiracSpinor&);
    /// Multiplication with scalar (dcomplex).
    const DiracSpinor operator * (const dcomplex&) const;
    /// Division by scalar.
    const DiracSpinor operator / (const dcomplex&) const;
    DiracSpinor& operator *= (const dcomplex&);
    DiracSpinor& operator /= (const dcomplex&);
    /// Friend operator for scalar*spinor.
    friend const DiracSpinor operator * (const dcomplex&, const DiracSpinor&);
    /// Matrix*spinor product.
    friend const DiracSpinor operator * (const DiracMatrix&, const DiracSpinor&);
    /// Scalar product.
    friend const dcomplex operator * (const AdDiracSpinor&, const DiracSpinor&);
    /// Diadic product.
    friend const DiracMatrix operator * (const DiracSpinor&, const AdDiracSpinor&);
    /// Equality check.
    bool operator == (const DiracSpinor&) const;
    /// Unequality check.
    bool operator != (const DiracSpinor&) const;
    /// Print the spinor to the given stream in the form '(x1,x2,x3,x4)'.
    friend std::ostream& operator << (std::ostream&, const DiracSpinor&);
    /// Read the spinor from the given stream. It must be in the form '(x1,x2,x3,x4)'.
    friend std::istream& operator >> (std::istream&, DiracSpinor&);
    /// Decide if any component is NaN.
    bool isnan() const;
    /// Dirac adjoint.
    const AdDiracSpinor adj() const;
    friend const AdDiracSpinor adj(const DiracSpinor&);

    friend class AdDiracSpinor;
    friend class DiracMatrix;
  private:
    /// Variables to store the spinor components.
    dcomplex x[4];
  };

  DiracSpinor U_(const Vectors::FourVector& p, halfint spin);
  DiracSpinor V_(const Vectors::FourVector& p, halfint spin);

  /**
     Represent an adjoint four-component (Dirac) spinor.
     Has the usual functionality (addition, multiplication with scalar, 
     inner product, adjoint, etc.) 
  */
  class AdDiracSpinor {
  public:
    /// Create a nullspinor.
    AdDiracSpinor();
    /// Copy constructor.
    AdDiracSpinor(const AdDiracSpinor&);
    /// Create a spinor with the given components.
    AdDiracSpinor(dcomplex,dcomplex,dcomplex,dcomplex);
    /// Create a spinor from two Pauli spinors.
    AdDiracSpinor(const AdPauliSpinor&, const AdPauliSpinor&);
    /// Destructor. Does nothing.
    ~AdDiracSpinor();
    /// Assignment operator.
    AdDiracSpinor& operator = (const AdDiracSpinor&);
    /// Positive sign.
    const AdDiracSpinor& operator + () const;
    /// Negative sign.
    const AdDiracSpinor operator - () const;
    /// Spinor addition.
    const AdDiracSpinor operator + (const AdDiracSpinor&) const;
    /// Spinor subtraction.
    const AdDiracSpinor operator - (const AdDiracSpinor&) const;
    /// Increment operator.
    AdDiracSpinor& operator += (const AdDiracSpinor&);
    /// Decrement operator.
    AdDiracSpinor& operator -= (const AdDiracSpinor&);
    /// Multiplication with scalar (dcomplex).
    const AdDiracSpinor operator * (const dcomplex&) const;
    /// Division by scalar.
    const AdDiracSpinor operator / (const dcomplex&) const;
    AdDiracSpinor& operator *= (const dcomplex&);
    AdDiracSpinor& operator /= (const dcomplex&);
    /// Friend operator for scalar*spinor.
    friend const AdDiracSpinor operator * (const dcomplex&, const AdDiracSpinor&);
    /// spinor*matrix product.
    friend const AdDiracSpinor operator * (const AdDiracSpinor&, const DiracMatrix&);
    /// Scalar product.
    friend const dcomplex operator * (const AdDiracSpinor&, const DiracSpinor&);
    /// Diadic product.
    friend const DiracMatrix operator * (const DiracSpinor&, const AdDiracSpinor&);
    /// Equality check.
    bool operator == (const AdDiracSpinor&) const;
    /// Unequality check.
    bool operator != (const AdDiracSpinor&) const;
    /// Print the spinor to the given stream in the form '(x1,x2,x3,x4)'.
    friend std::ostream& operator << (std::ostream&, const AdDiracSpinor&);
    /// Read the spinor from the given stream. It must be in the form '(x1,x2,x3,x4)'.
    friend std::istream& operator >> (std::istream&, AdDiracSpinor&);
    /// Decide if any component is NaN.
    bool isnan() const;
    /// Dirac adjoint.
    const DiracSpinor adj() const;
    friend const DiracSpinor adj(const AdDiracSpinor&);

    friend class DiracSpinor;
    friend class DiracMatrix;
  private:
    /// Variables to store the spinor components.
    dcomplex x[4];
  };

  AdDiracSpinor Ubar_(const Vectors::FourVector& p, halfint spin);
  AdDiracSpinor Vbar_(const Vectors::FourVector& p, halfint spin);

  /**
     Represent a Dirac matrix.
  */
  class DiracMatrix {
  public:
    /// Create a nullmatrix.
    DiracMatrix();
    /// Copy constructor.
    DiracMatrix(const DiracMatrix&);
    /// Create a Dirac matrix from components.
    DiracMatrix(dcomplex,dcomplex,dcomplex,dcomplex,
                dcomplex,dcomplex,dcomplex,dcomplex,
                dcomplex,dcomplex,dcomplex,dcomplex,
                dcomplex,dcomplex,dcomplex,dcomplex);
    /// Create a Dirac matrix from Pauli matrices.
    DiracMatrix(PauliMatrix,PauliMatrix,PauliMatrix,PauliMatrix);
    /// Destructor. Does nothing.
    ~DiracMatrix();
    /// Assignment operator.
    DiracMatrix& operator = (const DiracMatrix&);
    /// Adjoint
    const DiracMatrix adj() const;
    friend const DiracMatrix adj(const DiracMatrix&);
    /// Positive sign.
    const DiracMatrix& operator + () const;
    /// Negative sign.
    const DiracMatrix operator - () const;
    /// Matrix addition.
    const DiracMatrix operator + (const DiracMatrix&) const;
    /// Matrix subtraction.
    const DiracMatrix operator - (const DiracMatrix&) const;
    /// Increment operator.
    DiracMatrix& operator += (const DiracMatrix&);
    /// Decrement operator.
    DiracMatrix& operator -= (const DiracMatrix&);
    /// Matrix product.
    friend const DiracMatrix operator * (const DiracMatrix&, const DiracMatrix&);
    /// Matrix*spinor product.
    friend const DiracSpinor operator * (const DiracMatrix&, const DiracSpinor&);
    /// Spinor*matrix product.
    friend const AdDiracSpinor operator * (const AdDiracSpinor&, const DiracMatrix&);
    /// Diadic product.
    friend const DiracMatrix operator * (const DiracSpinor&, const AdDiracSpinor&);
    /// Multiplication with scalar (dcomplex).
    const DiracMatrix operator * (const dcomplex&) const;
    /// Division by scalar.
    const DiracMatrix operator / (const dcomplex&) const;
    DiracMatrix& operator *= (const dcomplex&);
    DiracMatrix& operator /= (const dcomplex&);
    /// Friend operator for scalar*matrix.
    friend const DiracMatrix operator * (const dcomplex&, const DiracMatrix&);
    /// Equality check.
    bool operator == (const DiracMatrix&) const;
    /// Unequality check.
    bool operator != (const DiracMatrix&) const;
    /// Print the matrix to the given stream in the form '((x11,x12,x13),(...),(...))'.
    friend std::ostream& operator << (std::ostream&, const DiracMatrix&);
    /// Print to a stream in the form of a table.
    void output(std::ostream&, int fieldWidth=20) const;
    /// Trace.
    dcomplex trace() const;
    friend dcomplex trace(const DiracMatrix&);
    friend bool isNaN(const DiracMatrix&);
    friend class DiracSpinor;
    friend class AdDiracSpinor;
  private:
    /// Variables to store the matrix components.
    dcomplex x[4][4];
    const DiracMatrix conj() const;
    const DiracMatrix transpose() const;
    const DiracMatrix trans_conj() const;
  };

  extern const DiracMatrix gamma0_;
  extern const DiracMatrix gamma1_;
  extern const DiracMatrix gamma2_;
  extern const DiracMatrix gamma3_;
  extern const DiracMatrix gamma5_;
  extern const DiracMatrix gamma_null;
  extern const DiracMatrix gamma_unit;

  class Gamma_mu {
  public:
    static const Gamma_mu& instance();
    const DiracMatrix& operator()(int) const;
    const DiracMatrix operator()(const Vectors::FourVector&) const;
  private:
    static const Gamma_mu gamma_;
    Gamma_mu();
    Gamma_mu(const Gamma_mu&);
  };

  extern const Gamma_mu& gamma_;

  class Sigma_mu_nu {
  public:
    static const Sigma_mu_nu& instance();
    const DiracMatrix operator()(int, int) const;
    const DiracMatrix operator()(int, const Vectors::FourVector&) const;
    const DiracMatrix operator()(const Vectors::FourVector&, int) const;
    const DiracMatrix operator()(const Vectors::FourVector&, const Vectors::FourVector&) const;
  private:
    static const Sigma_mu_nu sigma_;
    Sigma_mu_nu();
    Sigma_mu_nu(const Sigma_mu_nu&);
  };

  extern const Sigma_mu_nu& sigma_;

} // namespace feyntools

#endif
