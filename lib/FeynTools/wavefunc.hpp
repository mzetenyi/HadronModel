#ifndef WAVWFUNC_HPP
#define WAVWFUNC_HPP

#include <cassert>
#include <cstdlib>
#include "utils.hpp"
#include "Vectors.hpp"
#include "Spinors.hpp"
#include "MultiArray.hpp"

using namespace Vectors;
using namespace Spinors;

/**
   Represent the spin-1 polarization vector.
*/
class eps_ {
public:
  eps_(FourVector p);
  /**
     Returns the contravariant (upper index) component
     of the polarization vector.
     @param mu Lorentz index
     @param lambda polarization (2*spin z component)
  */
  dcomplex operator()(uint mu, halfint lambda);
private:
  MultiArray<dcomplex> val;

  /**
     Polarization threevector. This is the polarization vector for a spin-1 particle at rest.
     @param i index 
     @param lambda polarization (2*spin z coponent)
  */
  static dcomplex e_(int i, int lambda);
};

/**
   Polarization vector for bosons. Returns 1 for scalars (J=1).
   Returns the contravariant (upper index) component
   of the polarization vector for vector particles.
*/
class e_ {
public:
  /**
     Constructor.
     @param J total spin (J=0 and 1 are implemented)
     @param p four-momentum
  */
  e_(int J, FourVector p);
  /**
     Destructor. Delete dynamically allocated storage.
   */
  ~e_();
  /**
     Return components of the polarization vector.
     @param mu Lorentz index
     @param lambda polarization (spin z component)
  */
  const dcomplex operator()(int mu, int lambda);
  /**
     Returns the polarization fourvector as a FourVector object.
     @param lambda polarization (spin z component)

     \todo We need a complex FourVector
  */
  //const FourVector ...
  
private:
  const halfint J;
  eps_* eps_p;
};

class u_1h {
public:
  u_1h(FourVector p);
  //u_1h(Particle par);
  const DiracSpinor operator()(halfint spin);
private:
  MultiArray<DiracSpinor> val;
  void init(FourVector p);
};

class ubar_1h {
public:
  ubar_1h(FourVector p);
  //ubar_1h(Particle par);
  const AdDiracSpinor operator()(halfint spin);
private:
  MultiArray<AdDiracSpinor> val;
  void init(FourVector p);
};

class v_1h {
public:
  v_1h(FourVector p);
  //v_1h(Particle par);
  const DiracSpinor operator()(halfint spin);
private:
  MultiArray<DiracSpinor> val;
  void init(FourVector p);
};

/**
   The spin-3/2 spinor.
*/
class u_3h {
public:
  u_3h(FourVector p);
  //u_3h(Particle par);
  const DiracSpinor operator()(int mu, halfint spin);
private:
  MultiArray<DiracSpinor> val;
  void init(FourVector p);
};

/**
   The adjoint spin-3/2 spinor.
*/
class ubar_3h {
public:
  ubar_3h(FourVector p);
  //ubar_3h(Particle par);
  const AdDiracSpinor operator()(int mu, halfint spin);
private:
  MultiArray<AdDiracSpinor> val;
  void init(FourVector p);
};

/**
   The spin-5/2 spinor.
   @param mu Lorentz index
   @param nu Lorentz index
   @param p four-momentum
   @param s polarization (2*spin z coponent)
*/
class u_5h {
public:
  u_5h(FourVector p);
  //u_5h(Particle par);
  const DiracSpinor operator()(int mu, int nu, halfint spin);
private:
  MultiArray<DiracSpinor> val;
  void init(FourVector p);
};

/**
   The adjoint spin-5/2 spinor.
*/
class ubar_5h {
public:
  ubar_5h(FourVector p);
  //ubar_5h(Particle par);
  const AdDiracSpinor operator()(int mu, int nu, halfint spin);
private:
  MultiArray<AdDiracSpinor> val;
  void init(FourVector p);
};

/**
   Spinor amplitude for 1/2, 3/2 and 5/2 spin.

   \todo Do we need all constructors?
*/
class u_ {
public:
  /**
     Constructor.
     @param J total spin
     @param p four-momentum
  */
  u_(halfint J, FourVector p);
  /**
     Constructor.
     @param part a Particle object
  */
  //u_(Particle part);
  /**
     Destructor. Delete dynamically allocated storage.
  */
  ~u_();
  /** 
      Return the components of the spinor amplitude.
      @param mu "cumulated" Lorentz indices: mu = mu1 + 4*(mu2 + 4*(mu3 + ...) ).
      @param s polarization (spin z coponent)
   */
  const DiracSpinor operator()(int mu, halfint s);

  DiracSpinor operator()(halfint la, int mu1=0, int mu2=0);
private:
  const halfint J;
  u_1h* u_1hp;
  u_3h* u_3hp;
  u_5h* u_5hp;
};

/**
   Adjoint spinor amplitude for 1/2, 3/2 and 5/2 spin.
*/
class ubar_ {
public:
  /**
     Constructor.
     @param J total spin
     @param p four-momentum
  */
  ubar_(halfint J, FourVector p);
  /**
     Constructor.
     @param part a Particle object.
  */
  //ubar_(Particle part);
  /**
     Destructor. Delete dynamically allocated storage.
  */
  ~ubar_();
  /** 
      Return the components of the spinor amplitude.
      @param mu "cumulated" Lorentz indices: mu = mu1 + 4*(mu2 + 4*(mu3 + ...) ).
      @param s polarization (spin z coponent)
   */
  const AdDiracSpinor operator()(int mu, halfint s);
  AdDiracSpinor operator()(halfint la, int mu1=0, int mu2=0);
private:
  const halfint J;
  ubar_1h* ubar_1hp;
  ubar_3h* ubar_3hp;
  ubar_5h* ubar_5hp;
};

/**
  The spin-1 projector P(mu,nu) = Sum_lambda e_(mu,p,lambda)*e_(nu,p,lambda) = g_(mu,nu) - p(mu)*p(nu)/m^2.
*/
class pro1 {
public:
  //  pro1(Particle par);
  const double operator()(int mu, int nu);
private:
  MultiArray<double> val;
};

/**
   The spin-3/2 projector. Indices mu and nu are upper.
*/
class pro3h {
public:
  //  pro3h(Particle par);
  const DiracMatrix operator()(int mu, int nu);
private:
  MultiArray<DiracMatrix> val;
};


/**
   The spin-5/2 projector. Indices mu1, mu2, nu1 and nu2 are upper.
*/
class pro5h {
public:
  //  pro5h(Particle par);
  const DiracMatrix operator()(int mu1, int mu2, int nu1, int nu2);
private:
  MultiArray<DiracMatrix> val;
  static double _G(double m, int mu, int nu, FourVector p);
  static DiracMatrix _T(double m, int mu, int nu, FourVector p);
};

#endif // WAVWFUNC_HPP
