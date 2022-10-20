#include "wavefunc.hpp"

#include <cassert>
#include <cstdlib>

using namespace std;

eps_::eps_(FourVector p) : val(idx_lor, idx_s1) {
  double m2 = p * p;
  if (m2 < 0) {
    if (-m2 / (p(0) * p(0)) > 0.000001) {
      cerr << "spacelike fourmomentum in eps_: " << p << endl;
      exit(0);
    }
    m2 = 0;
  }
  if (m2 / (p(0) * p(0)) <
      0.000001) {  // zero mass particle -- in radiation gauge, i.e. e_(0)=0
    const ThreeVector pp = p.spacial();
    const ThreeVector ex(1, 0, 0);
    const ThreeVector ey(0, 1, 0);
    ThreeVector e1 =
        ((ex * pp < sqrt(pp * pp / 2)) ? crossprod(pp, ex) : crossprod(pp, ey));
    e1 /= sqrt(e1 * e1);
    ThreeVector e2 = crossprod(pp, e1);
    e2 /= sqrt(e2 * e2);
    const FourVector eps0(1, 0, 0, 0);
    const FourVector eps1(0, e1);
    const FourVector eps2(0, e2);
    // const FourVector eps3 = (p - (p*eps0)*eps0)/sqrt((p*eps0)*(p*eps0) -
    // p*p);
    for (halfint mu(_0); mu < _4; mu++) {
      val(mu, _1) = -(eps1(mu) + i_ * eps2(mu)) / sqrt(2.);
      val(mu, _0) = 0;
      val(mu, -_1) = (eps1(mu) - i_ * eps2(mu)) / sqrt(2.);
    }
  } else {
    // massive particle:
    double m = sqrt(m2);
    for (halfint lambda(-_1); lambda <= _1; lambda++) {
      dcomplex ep(0);
      for (int i(1); i <= 3; i++) {
        ep += e_(i, lambda) * p(i);
      }
      val(_0, lambda) = ep / m;
      for (halfint i(_1); i <= _3; i++) {
        val(i, lambda) = e_(i, lambda) + ep * p(i) / (m * (p(0) + m));
      }
    }
  }
}

/**
   Returns the contravariant (upper index) component
   of the polarization vector.
   @param mu Lorentz index
   @param lambda polarization (2*spin z component)
*/
dcomplex eps_::operator()(uint mu, halfint lambda) {
  assert(isIndex4(mu));
  assert(is_spin(_1, lambda));
  return val(_1 * mu, lambda);
}

/**
   Polarization threevector. This is the polarization vector for a spin-1
   particle at rest.
   @param i index
   @param lambda polarization (2*spin z coponent)
*/
dcomplex eps_::e_(int i, int lambda) {
  assert(isIndex3(i));
  assert(lambda == -1 or lambda == 0 or lambda == 1);
  if (lambda == 0) {
    if (i == 3) {
      return 1;
    }
    return 0;
  }
  if (i == 3) {
    return 0;
  }
  if (i == 2) {
    return -i_ / sqrt(2.);
  }
  if (lambda == 1) {
    return -1. / sqrt(2.);
  }
  return 1. / sqrt(2.);
}

/**
   Constructor.
   @param J total spin (J=0 and 1 are implemented)
   @param p four-momentum
*/
e_::e_(int J, FourVector p) : J(J), eps_p(NULL) {
  assert((J == 0) or (J == 1));
  if (J == 1) {
    eps_p = new eps_(p);
  }
}
/**
     Destructor. Delete dynamically allocated storage.
*/
e_::~e_() {
  if (eps_p) {
    delete eps_p;
  }
}

/**
   Return components of the polarization vector.
   @param mu Lorentz index
   @param lambda polarization (spin z component)
*/
const dcomplex e_::operator()(int mu, int lambda) {
  assert(is_spin(J, _1 * lambda));
  if (J == 0) {
    assert(mu == 0);
    return 1;
  }
  if (J == 1) {
    assert(0 <= mu and mu < 4);
    return (*eps_p)(mu, _1 * lambda);
  }
  return 0;
}

u_1h::u_1h(FourVector p) : val(idx_s1h) { init(p); }

/*
u_1h::u_1h(Particle par) : val(idx_s1h) {
  init(par.fourmom());
}
*/

void u_1h::init(FourVector p) {
  double m2 = p * p;
  if (m2 < 0) {
    if (-m2 / (p(0) * p(0)) > 0.000001) {
      cerr << "spacelike fourmomentum in u_1h: " << p << endl;
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
  for (halfint spin(-half); spin <= half; spin++) {
    PauliSpinor khi = Khi(spin);
    val(spin) = DiracSpinor(norm * khi, 1. / norm * (tau_p * khi));
  }
}

const DiracSpinor u_1h::operator()(halfint spin) {
  assert(is_spin(half, spin));
  return val(spin);
}

ubar_1h::ubar_1h(FourVector p) : val(idx_s1h) { init(p); }

/*
ubar_1h::ubar_1h(Particle par) : val(idx_s1h) {
  init(par.fourmom());
}
*/

void ubar_1h::init(FourVector p) {
  u_1h u(p);
  for (halfint spin(-half); spin <= half; spin++) {
    val(spin) = adj(u(spin));
  }
}

const AdDiracSpinor ubar_1h::operator()(halfint spin) {
  assert(is_spin(half, spin));
  return val(spin);
}

v_1h::v_1h(FourVector p) : val(idx_s1h) { init(p); }

/*
v_1h::v_1h(Particle par) : val(idx_s1h) {
  init(par.fourmom());
}
*/

void v_1h::init(FourVector p) {
  assert(p(0) > 0);
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
  for (halfint spin(-half); spin <= half; spin++) {
    PauliSpinor khi = Khi(-spin);
    val(spin) = DiracSpinor(spin / norm * (tau_p * khi), spin * norm * khi);
  }
}

const DiracSpinor v_1h::operator()(halfint spin) {
  assert(is_spin(half, spin));
  return val(spin);
}

u_3h::u_3h(FourVector p) : val(idx_lor, idx_s3h) { init(p); }

/*
u_3h::u_3h(Particle par) : val(idx_lor,idx_s3h) {
  init(par.fourmom());
}
*/

void u_3h::init(FourVector p) {
  double m2 = p * p;
  if (m2 < 0) {
    if (-m2 / (p(0) * p(0)) > 0.000001) {
      cerr << "spacelike fourmomentum in u_(int mu, FourVector p, int s): " << p
           << endl;
      exit(0);
    }
    m2 = 0;
  }
  if (m2 / (p(0) * p(0)) < 0.000001) {  // zero mass particle
    cerr << "u_(int mu, FourVector p, int s): zero mass particle - not "
            "implemented."
         << endl;
    exit(0);
  }
  u_1h u_(p);
  eps_ e_(p);
  for (halfint mu(_0); mu < _4; mu++) {
    val(mu, 3 * half) = e_(mu, _1) * u_(half);
    val(mu, half) = 1. / sqrt(3.) * e_(mu, _1) * u_(-half) +
                    sqrt(2. / 3.) * e_(mu, _0) * u_(half);
    val(mu, -half) = sqrt(2. / 3.) * e_(mu, _0) * u_(-half) +
                     1. / sqrt(3.) * e_(mu, -_1) * u_(half);
    val(mu, -3 * half) = e_(mu, -_1) * u_(-half);
  }
}

const DiracSpinor u_3h::operator()(int mu, halfint spin) {
  assert(isIndex4(mu));
  assert(is_spin(3 * half, spin));
  return val(_1 * mu, spin);
}

ubar_3h::ubar_3h(FourVector p) : val(idx_lor, idx_s3h) { init(p); }

/*
ubar_3h::ubar_3h(Particle par) : val(idx_lor,idx_s3h) {
  init(par.fourmom());
}
*/

void ubar_3h::init(FourVector p) {
  u_3h u(p);
  for (halfint mu(_0); mu < _4; mu++) {
    for (halfint spin(-3 * half); spin <= 3 * half; spin++) {
      val(mu, spin) = adj(u(mu, spin));
    }
  }
}

const AdDiracSpinor ubar_3h::operator()(int mu, halfint spin) {
  assert(isIndex4(_1 * mu));
  assert(is_spin(3 * half, spin));
  return val(_1 * mu, spin);
}

u_5h::u_5h(FourVector p) : val(idx_lor, idx_lor, idx_s5h) { init(p); }

/*
u_5h::u_5h(Particle par) : val(idx_lor,idx_lor,idx_s5h) {
  init(par.fourmom());
}
*/

void u_5h::init(FourVector p) {
  double m2 = p * p;
  if (m2 < 0) {
    if (-m2 / (p(0) * p(0)) > 0.000001) {
      cerr << "spacelike fourmomentum in u_(int mu, int nu, FourVector p, int "
              "s): "
           << p << endl;
      exit(0);
    }
    m2 = 0;
  }
  if (m2 / (p(0) * p(0)) < 0.000001) {  // zero mass particle
    cerr << "u_(int mu, int nu, FourVector p, int s): zero mass particle - not "
            "implemented."
         << endl;
    exit(0);
  }
  u_3h u_(p);
  eps_ e_(p);
  for (halfint mu(_0); mu < _4; mu++) {
    for (halfint nu(_0); nu < _4; nu++) {
      val(mu, nu, 5 * half) = e_(mu, _1) * u_(nu, 3 * half);
      val(mu, nu, 3 * half) = sqrt(2. / 5.) * e_(mu, _0) * u_(nu, 3 * half) +
                              sqrt(3. / 5.) * e_(mu, _1) * u_(nu, half);
      val(mu, nu, 1 * half) = 1. / sqrt(10.) * e_(mu, -_1) * u_(nu, 3 * half) +
                              sqrt(3. / 5.) * e_(mu, _0) * u_(nu, half) +
                              sqrt(3. / 10.) * e_(mu, _1) * u_(nu, -half);
      val(mu, nu, -1 * half) = sqrt(3. / 10.) * e_(mu, -_1) * u_(nu, half) +
                               sqrt(3. / 5.) * e_(mu, _0) * u_(nu, -half) +
                               1. / sqrt(10.) * e_(mu, _1) * u_(nu, -3 * half);
      val(mu, nu, -3 * half) = sqrt(3. / 5.) * e_(mu, -_1) * u_(nu, -half) +
                               sqrt(2. / 5.) * e_(mu, _0) * u_(nu, -3 * half);
      val(mu, nu, -5 * half) = e_(mu, -_1) * u_(nu, -3 * half);
    }
  }
}

const DiracSpinor u_5h::operator()(int mu, int nu, halfint spin) {
  assert(isIndex4(mu));
  assert(isIndex4(nu));
  assert(is_spin(5 * half, spin));
  return val(_1 * mu, _1 * nu, spin);
}

ubar_5h::ubar_5h(FourVector p) : val(idx_lor, idx_lor, idx_s5h) { init(p); }

/*
ubar_5h::ubar_5h(Particle par) : val(idx_lor,idx_lor,idx_s5h) {
  init(par.fourmom());
}
*/

void ubar_5h::init(FourVector p) {
  cerr << "ubar_5h::init(FourVector p) is not correctly implemented!!!" << endl;
  u_3h u(p);
  for (halfint mu(_0); mu < _4; mu++) {
    for (halfint nu(_0); nu < _4; nu++) {
      for (halfint spin(-5 * half); spin <= 5 * half; spin++) {
        val(mu, spin) = adj(u(mu, spin));
      }
    }
  }
}

const AdDiracSpinor ubar_5h::operator()(int mu, int nu, halfint spin) {
  assert(isIndex4(mu));
  assert(isIndex4(nu));
  assert(is_spin(5 * half, spin));
  return val(_1 * mu, _1 * nu, spin);
}

/**
   Constructor.
   @param J total spin
   @param p four-momentum
*/
u_::u_(halfint J, FourVector p) : J(J), u_1hp(NULL), u_3hp(NULL), u_5hp(NULL) {
  assert(p.future());
  if (J == half) {
    u_1hp = new u_1h(p);
  } else if (J == 3 * half) {
    u_3hp = new u_3h(p);
  } else if (J == 5 * half) {
    u_5hp = new u_5h(p);
  } else {
    cerr << "u_(J,p):  J = " << J << " not implemented." << endl;
    exit(0);
  }
}
/**
   Constructor.
   @param part a Particle object

   \todo 0 mass particles are excluded by assertions.
*/
/*
u_::u_(Particle part) : J(part.spin()), u_1hp(NULL), u_3hp(NULL), u_5hp(NULL) {
  assert(part.fourmom().future());
  if (J==half) {
    u_1hp = new u_1h(part.fourmom());
  } else if (J==3*half) {
    u_3hp = new u_3h(part.fourmom());
  } else if (J==5*half) {
    u_5hp = new u_5h(part.fourmom());
  } else {
    cerr << "u_(Particle):  invalid or unimplemented spin of baryon, J = " << J
<< endl; exit(0);
  }
}
*/

/**
     Destructor. Delete dynamically allocated storage.
*/
u_::~u_() {
  if (u_1hp) {
    delete u_1hp;
  }
  if (u_3hp) {
    delete u_3hp;
  }
  if (u_5hp) {
    delete u_5hp;
  }
}
/**
    Return the components of the spinor amplitude.
    @param mu "cumulated" Lorentz indices: mu = mu1 + 4*(mu2 + 4*(mu3 + ...) ).
    @param s polarization (spin z coponent)
*/
const DiracSpinor u_::operator()(int mu, halfint s) {
  assert(mu >= 0);
  if (J == half) {
    assert(mu == 0);
    return (*u_1hp)(s);
  } else if (J == 3 * half) {
    assert(mu < 4);
    return (*u_3hp)(mu, s);
  } else if (J == 5 * half) {
    assert(mu < 16);
    int mu1 = mu % 4;
    int mu2 = mu / 4;
    return (*u_5hp)(mu1, mu2, s);
  } else {
    cerr << "u_(mu,s):  J = " << J << " not implemented." << endl;
    exit(0);
  }
  DiracSpinor ret;
  return ret;  // never get here
}

DiracSpinor u_::operator()(halfint la, int mu1, int mu2) {
  return operator()(mu1 + 4 * mu2, la);
}

/**
   Constructor.
   @param J total spin
   @param p four-momentum
*/
ubar_::ubar_(halfint J, FourVector p)
    : J(J), ubar_1hp(NULL), ubar_3hp(NULL), ubar_5hp(NULL) {
  assert(p.future());
  if (J == half) {
    ubar_1hp = new ubar_1h(p);
  } else if (J == 3 * half) {
    ubar_3hp = new ubar_3h(p);
  } else if (J == 5 * half) {
    ubar_5hp = new ubar_5h(p);
  } else {
    cerr << "ubar_(J,p):  J = " << J << " not implemented." << endl;
    exit(0);
  }
}
/**
   Constructor.
   @param part a Particle object.
*/
/*
ubar_::ubar_(Particle part) : J(part.spin()), ubar_1hp(NULL), ubar_3hp(NULL),
ubar_5hp(NULL) { assert(part.fourmom().future()); if (J==half) { ubar_1hp = new
ubar_1h(part.fourmom()); } else if (J==3*half) { ubar_3hp = new
ubar_3h(part.fourmom()); } else if (J==5*half) { ubar_5hp = new
ubar_5h(part.fourmom()); } else { cerr << "ubar_(Particle):  invalid or
unimplemented spin of baryon, J = " << J << endl; exit(0);
  }
}
*/

/**
     Destructor. Delete dynamically allocated storage.
*/
ubar_::~ubar_() {
  if (ubar_1hp) {
    delete ubar_1hp;
  }
  if (ubar_3hp) {
    delete ubar_3hp;
  }
  if (ubar_5hp) {
    delete ubar_5hp;
  }
}
/**
    Return the components of the spinor amplitude.
    @param mu "cumulated" Lorentz indices: mu = mu1 + 4*(mu2 + 4*(mu3 + ...) ).
    @param s polarization (spin z coponent)
*/
const AdDiracSpinor ubar_::operator()(int mu, halfint s) {
  assert(mu >= 0);
  if (J == half) {
    assert(mu == 0);
    return (*ubar_1hp)(s);
  } else if (J == 3 * half) {
    assert(mu < 4);
    return (*ubar_3hp)(mu, s);
  } else if (J == 5 * half) {
    assert(mu < 16);
    int mu1 = mu % 4;
    int mu2 = mu / 4;
    return (*ubar_5hp)(mu1, mu2, s);
  } else {
    cerr << "ubar_(mu,s):  J = " << J << " not implemented." << endl;
    exit(0);
  }
  AdDiracSpinor ret;
  return ret;  // never get here
}

AdDiracSpinor ubar_::operator()(halfint la, int mu1, int mu2) {
  return operator()(mu1 + 4 * mu2, la);
}

/*
pro1::pro1(Particle par) : val(idx_lor,idx_lor) {
  FourVector p = par.fourmom();
  double m = par.mass();
  for (int mu(0); mu<4; mu++) {
    for (int nu(0); nu<4; nu++) {
      val(mu,nu) =  - g_(mu,nu) + p(mu)*p(nu)/(m*m);
    }
  }
}
*/

const double pro1::operator()(int mu, int nu) {
  assert(idx_lor.match(_1 * mu));
  assert(idx_lor.match(_1 * nu));
  return val(_1 * mu, _1 * nu);
}

/*
pro3h::pro3h(Particle par) : val(idx_lor,idx_lor) {
  FourVector p = par.fourmom();
  double m = par.mass();
  for (int mu(0); mu<4; mu++) {
    for (int nu(0); nu<4; nu++) {
      val(mu,nu) +=
        - g_(mu,nu) * gamma_unit
        + gamma_(mu)*gamma_(nu)/3.
        + 2.*p(mu)*p(nu)/(3.*m*m) * gamma_unit
        - (p(mu)*gamma_(nu) - p(nu)*gamma_(mu))/(3.*m);
    }
  }
}
*/

const DiracMatrix pro3h::operator()(int mu, int nu) {
  assert(idx_lor.match(_1 * mu));
  assert(idx_lor.match(_1 * nu));
  return val(_1 * mu, _1 * nu);
}

/*
pro5h::pro5h(Particle par) : val(idx_lor,idx_lor,idx_lor,idx_lor) {
  FourVector p = par.fourmom();
  double m = par.mass();
  MultiArray<double> G(idx_lor,idx_lor);
  MultiArray<DiracMatrix> T(idx_lor,idx_lor);
  for (int mu(0); mu<4; mu++) {
    for (int nu(0); nu<4; nu++) {
      G(mu,nu) = _G(m,mu,nu,p);
      T(mu,nu) = _T(m,mu,nu,p);
    }
  }
  for (int mu1(0); mu1<4; mu1++) {
    for (int mu2(0); mu2<4; mu2++) {
      for (int nu1(0); nu1<4; nu1++) {
        for (int nu2(0); nu2<4; nu2++) {
          val(mu1,mu2,nu1,nu2) +=
            - 3./10. * (G(mu1,nu1)*G(mu2,nu2) + G(mu1,nu2)*G(mu2,nu1)) *
gamma_unit
            + 1./5. * G(mu1,mu2)*G(nu1,nu2) * gamma_unit
            + 1./10. * (T(mu1,nu1)*G(mu2,nu2) + T(mu2,nu2)*G(mu1,nu1) +
T(mu1,nu2)*G(mu2,nu1) + T(mu2,nu1)*G(mu1,nu2));
        }
      }
    }
  }
}
*/

const DiracMatrix pro5h::operator()(int mu1, int mu2, int nu1, int nu2) {
  assert(idx_lor.match(_1 * mu1));
  assert(idx_lor.match(_1 * mu2));
  assert(idx_lor.match(_1 * nu1));
  assert(idx_lor.match(_1 * nu2));
  return val(_1 * mu1, _1 * mu2, _1 * nu1, _1 * nu2);
}

double pro5h::_G(double m, int mu, int nu, FourVector p) {
  assert(idx_lor.match(_1 * mu));
  assert(idx_lor.match(_1 * nu));
  return -g_(mu, nu) + p(mu) * p(nu) / (m * m);
}

DiracMatrix pro5h::_T(double m, int mu, int nu, FourVector p) {
  return -1. / 2. * (gamma_(mu) * gamma_(nu) - gamma_(nu) * gamma_(mu)) +
         1. / (2. * m) * p(mu) *
             (gamma_(p) * gamma_(nu) - gamma_(nu) * gamma_(p)) -
         1. / (2. * m) * p(nu) *
             (gamma_(p) * gamma_(mu) - gamma_(mu) * gamma_(p));
}
