#include "FeynTools.hpp"

using namespace units_GeV;

DiracMatrix pro12analytic(FourVector p, double m) {
  double minv = sqrt(p * p);
  return gamma_(p) + minv * gamma_unit;
}

DiracMatrix pro12numeric(FourVector p, double m) {
  u_ u(half, p);
  ubar_ ubar(half, p);
  DiracMatrix ret = gamma_null;
  for (halfint la : {-half, half}) {
    ret += u(la) * ubar(la);
  }
  return ret;
}

DiracMatrix pro32analytic(FourVector p, double m, int mu, int nu,
                          int version = 0) {
  assert(isIndex4(mu));
  assert(isIndex4(nu));
  double m2 = m * m;
  double minv2 = p * p;
  double minv = sqrt(minv2);
  if (version == 0) {
    return (gamma_(p) + m * gamma_unit) *
           ((-g_(mu, nu) + 2. / (3. * m2) * p(mu) * p(nu)) * gamma_unit +
            gamma_(mu) * gamma_(nu) / 3. -
            (gamma_(mu) * p(nu) - gamma_(nu) * p(mu)) / (3. * m));
  }
  return (gamma_(p) + minv * gamma_unit) *
         ((-g_(mu, nu) + 2. / (3. * minv2) * p(mu) * p(nu)) * gamma_unit +
          gamma_(mu) * gamma_(nu) / 3. -
          (gamma_(mu) * p(nu) - gamma_(nu) * p(mu)) / (3. * minv));
}

DiracMatrix pro32numeric(FourVector p, double m, int mu, int nu) {
  assert(isIndex4(mu));
  assert(isIndex4(nu));
  double m2 = m * m;
  u_ u(3 * half, p);
  ubar_ ubar(3 * half, p);
  DiracMatrix ret(gamma_null);
  for (halfint la : {-3 * half, -half, half, 3 * half}) {
    ret += u(la, mu) * ubar(la, nu);
  }
  return ret;
}

int main(int argc, char** argv) {
  Config::load(argc, argv);
  double m = getParam<double>("m", 1.5 * GeV);
  double minv = getParam<double>("minv", 1.5 * GeV);
  double pabs = getParam<double>("pabs", 0.);
  double costh = getParam<double>("costh", 1.);
  double sinth = sqrt(1. - costh * costh);
  double phi = getParam<double>("phi", 0.);
  double E = sqrt(minv * minv + pabs * pabs);
  FourVector p(E, pabs * costh, pabs * sinth * cos(phi),
               pabs * sinth * sin(phi));
  for (int mu : {0, 1, 2, 3}) {
    for (int nu : {0, 1, 2, 3}) {
      cout << "mu = " << mu << "   nu = " << nu << endl;
      cout << setw(4) << "i" << setw(4) << "j" << setw(18) << "analytic0"
           << setw(18) << "analytic1" << setw(18) << "numeric" << setw(18)
           << "diff" << endl;
      DiracMatrix pro32anal0 = pro32analytic(p, m, mu, nu, 0);
      DiracMatrix pro32anal1 = pro32analytic(p, m, mu, nu, 1);
      DiracMatrix pro32num = pro32numeric(p, m, mu, nu);
      for (int i : {0, 1, 2, 3}) {
        for (int j : {0, 1, 2, 3}) {
          cout << setw(4) << i << setw(4) << j << setw(18) << pro32anal0(i, j)
               << setw(18) << pro32anal1(i, j) << setw(18) << pro32num(i, j)
               << setw(18) << pro32num(i, j) - pro32anal1(i, j) << endl;
        }
      }
      cout << endl;
    }
  }
  cout << endl << "1/2 projector:" << endl;
  DiracMatrix pro12anal = pro12analytic(p, m);
  DiracMatrix pro12num = pro12numeric(p, m);
  cout << setw(4) << "i" << setw(4) << "j" << setw(18) << "analytic" << setw(18)
       << "numeric" << endl;
  for (int i : {0, 1, 2, 3}) {
    for (int j : {0, 1, 2, 3}) {
      cout << setw(4) << i << setw(4) << j << setw(18) << pro12anal(i, j)
           << setw(18) << pro12num(i, j) << endl;
    }
  }
  return 0;
}