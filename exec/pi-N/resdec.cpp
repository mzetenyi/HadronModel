#include "resdec.hpp"

//#include <iostream>
//#include <iomanip>

#include "utils.hpp"
#include "Vectors.hpp"
using namespace Vectors;
#include "Spinors.hpp"
using namespace Spinors;
#include "wavefunc.hpp"
#include "halfint.hpp"
#include "Config.hpp"
#include "Vrancx.hpp"
#include "Gamma.hpp"
#include "cutoffs.hpp"

using namespace std;

double Gsigma_pipi(double g, double mm) {
  const double msi = Config::get<double>("sigma.mass");
  const double mpi = Config::get<double>("pi_pm.mass");
  double pf = sqrt(mm*mm/4. - mpi*mpi);
  return 9.*g*g/(32.*pi_*mpi*mpi) * pf/(mm*msi) * pow(mm*mm/2. - mpi*mpi,2);
}

double GammaR_Npi(string res, double mr) {
  const double g = Config::get<double>(res+".g0");
  const halfint spin = Config::get<halfint>(res+".spin");
  const int parity = Config::get<int>(res+".parity");
  const halfint spinParity = spin*parity;
  const halfint isospin = Config::get<halfint>(res+".isospin");
  const double mn = Config::get<double>("Nucleon.mass");
  const double mpi_pm = Config::get<double>("pi_pm.mass");
  const double mpi_0 = Config::get<double>("pi_0.mass");
  //const double mpi = (2.*mpi_pm + mpi_0)/3.;
  const double mpi = mpi_pm;

  double npol = 2.*habs(spin) + 1.;
  double isofac = ((isospin==half) ? 3. : 1.);

  double pabs = sqrt(lambda(mr*mr,mn*mn,mpi*mpi))/(2.*mr);
  double EN = (mr*mr + mn*mn - mpi*mpi)/(2.*mr);
  double Epi = (mr*mr - mn*mn + mpi*mpi)/(2.*mr);
  FourVector pR(mr,0,0,0);
  FourVector pN(EN,0,0,-pabs);
  FourVector q(Epi,0,0,pabs);

  ubar_1h ubar(pN);

  double MSQR(0);

  //cerr << res << " -> Npi:  " << endl;

  // cutoff:
  double cutoff(1);
  if (Config::exists("cutoff")) {
    cutoff = cutoff_RNpi(res, mr);
  }
  
  if (spin==half) {
    u_1h u(pR);
    for (halfint laN : {-half,half}) {
      for (halfint laR : {-half,half}) {
        dcomplex M = cutoff * ubar(laN) * vertex1hNpi(g,spinParity,-q) * u(laR);
        //cerr << "  M = " << M << endl;
        MSQR += real(M*conj(M));
      }
    }
  } else if (spin==3*half) {
    u_3h u(pR);
    for (halfint laN : {-half,half}) {
      for (halfint laR(-3*half); laR<=3*half; laR+=_1) {
        dcomplex M(0);
        for (uint mu : {0,1,2,3}) {
          M += cutoff * ubar(laN) * vertex3hNpi(g,spinParity,mu,pR,-q) * sign_(mu) * u(mu,laR);
        }
        //cerr << "  M = " << M << endl;
        MSQR += real(M*conj(M));
      }
    }
  } else if (spin==5*half) {
    u_5h u(pR);
    for (halfint laN : {-half,half}) {
      for (halfint laR(-5*half); laR<=5*half; laR+=_1) {
        dcomplex M(0);
        for (uint mu : {0,1,2,3}) {
          for (uint nu : {0,1,2,3}) {
            M += cutoff * ubar(laN) * vertex5hNpi(g,spinParity,mu,nu,pR,-q) * sign_(mu)*sign_(nu) * u(mu,nu,laR);
          }
        }
        //cerr << "  M = " << M << endl;
        MSQR += real(M*conj(M));
      }
    }
  }

  //cerr << "  MSQR = " << MSQR << endl;

  return 1./(8.*pi_) * pabs/(mr*mr) * isofac/npol * MSQR;
}

double dGR_Nrho_dm(string res, double mr, double mm) { // mm: mass of rho
  const double g1 = Config::get<double>(res+".g1");
  const double g2(0);
  const double g3(0);
  const halfint spin = Config::get<halfint>(res+".spin");
  const int parity = Config::get<int>(res+".parity");
  const halfint spinParity = spin*parity;
  const halfint isospin = Config::get<halfint>(res+".isospin");
  const double mn = Config::get<double>("Nucleon.mass");
  const double mrho = Config::get<double>("rho.mass");

  double npol = 2.*habs(spin) + 1.;
  double isofac = ((isospin==half) ? 3. : 1.);

  double pabs = sqrt(lambda(mr*mr,mn*mn,mm*mm))/(2.*mr);
  double EN = (mr*mr + mn*mn - mm*mm)/(2.*mr);
  double Erho = (mr*mr - mn*mn + mm*mm)/(2.*mr);
  FourVector pR(mr,0,0,0);
  FourVector pN(EN,0,0,-pabs);
  FourVector k(Erho,0,0,pabs);

  ubar_1h ubar(pN);
  eps_ eps(k);

  // Rapp cutoff:
  double cutoff(1);
  // cutoff turned off for R N rho vertex:
  // if (Config::exists("cutoffRapp")) {
  //   cutoff = cutoffRapp(spin, pabs*pabs);
  // }
  
  double MSQR(0);

  if (spin==half) {
    u_1h u(pR);
    for (halfint laN : {-half,half}) {
      for (halfint larho : {-_1,_0,_1}) {
        for (halfint laR : {-half,half}) {
          dcomplex M(0);
          for (uint si : {0,1,2,3}) {
            M += cutoff * ubar(laN) * vertex1hNrho(g1,spinParity,pR,si,-k) * u(laR) * eps(si,larho) * sign_(si);
          }
          MSQR += real(M*conj(M));
        }
      }
    }
  } else if (spin==3*half) {
    u_3h u(pR);
    for (halfint laN : {-half,half}) {
      for (halfint larho : {-_1,_0,_1}) {
        for (halfint laR(-3*half); laR<=3*half; laR+=_1) {
          dcomplex M(0);
          for (uint mu : {0,1,2,3}) {
            for (uint si : {0,1,2,3}) {
              M += cutoff * ubar(laN) * vertex3hNrho(g1,g2,g3,spinParity,mu,pR,si,-k) * sign_(mu)*sign_(si) * u(mu,laR) * eps(si,larho);
            }
          }
          MSQR += real(M*conj(M));
        }
      }
    }
  } else if (spin==5*half) {
    u_5h u(pR);
    for (halfint laN : {-half,half}) {
      for (halfint larho : {-_1,_0,_1}) {
        for (halfint laR(-5*half); laR<=5*half; laR+=_1) {
          dcomplex M(0);
          for (uint mu : {0,1,2,3}) {
            for (uint nu : {0,1,2,3}) {
              for (uint si : {0,1,2,3}) {
                M += cutoff * ubar(laN) * vertex5hNrho(g1,g2,g3,spinParity,mu,nu,pR,si,-k) * sign_(mu)*sign_(nu)*sign_(si) * u(mu,nu,laR) * eps(si,larho);
              }
            }
          }
          MSQR += real(M*conj(M));
        }
      }
    }
  }

  double Gamma = Gamma_rho(mm);

  //return 1./(4.*pi_*pi_) * pabs/(mr*mr) * isofac/npol * MSQR * (mrho*mrho * Gamma)/(POW<2>(mrho*mrho-mm*mm) + mm*mm*Gamma*Gamma);
  return 1./(4.*pi_*pi_) * pabs/(mr*mr) * isofac/npol * MSQR * (mm*mrho * Gamma)/(POW<2>(mrho*mrho-mm*mm) + mm*mm*Gamma*Gamma);
  //return 1./(4.*pi_*pi_) * pabs/(mr*mr) * isofac/npol * MSQR * (mm*mrho * Gamma)/(POW<2>(mrho*mrho-mm*mm) + mrho*mrho*Gamma*Gamma); // verion used for publ.
}

double GammaR_Nrho(string res, double mr) {
  const double mn = Config::get<double>("Nucleon.mass");
  const double mpi = Config::get<double>("pi_pm.mass");
  const double mmin(2.*mpi);
  const double mmax(mr-mn);
  const int N(100);
  const double dm = (mmax-mmin)/N;
  double sum(0);
  for (double m = mmin+dm/2.; m<mmax; m+=dm) {
    sum += dGR_Nrho_dm(res, mr, m);
  }
  return sum*dm;
}

double GammaR_Ngamma(string res, double mr, string VMDversion) { // VMD1 is Kroll-Lee-Zumino, VMD2 is the other
  double g1; // RNrho coupling for VMD1, RNgamma coupling otherwise
  bool VMD2(VMDversion=="VMD2");
  if (VMD2) {
    g1 = Config::get<double>(res+".gVMD2");
  } else {
    g1 = Config::get<double>(res+".gngamma");
  }
  const double g2(0);
  const double g3(0);
  const double grho = Config::get<double>("grho");
  const halfint spin = Config::get<halfint>(res+".spin");
  const int parity = Config::get<int>(res+".parity");
  const halfint spinParity = spin*parity;
  const halfint isospin = Config::get<halfint>(res+".isospin");
  const double mn = Config::get<double>("Nucleon.mass");
  const double mrho = Config::get<double>("rho.mass");

  double npol = 2.*habs(spin) + 1.;
  double isofac = 1.;

  //  double pabs = sqrt(lambda(mr*mr,mn*mn,0))/(2.*mr);
  double pabs = (mr*mr-mn*mn)/(2.*mr);
  double EN = (mr*mr + mn*mn)/(2.*mr);
  double Egamma = pabs;
  FourVector pR(mr,0,0,0);
  FourVector pN(EN,0,0,-pabs);
  FourVector k(Egamma,0,0,pabs);

  // Rapp cutoff:
  double cutoff(1);
  if (Config::exists("cutoffRapp")) {
    cutoff = cutoffRapp(spin, pabs*pabs);
  }

  ubar_1h ubar(pN);
  eps_ eps(k);

  double MSQR(0);

  if (spin==half) {
    u_1h u(pR);
    for (halfint laN : {-half,half}) {
      for (halfint lag : {-_1,_0,_1}) {
        for (halfint laR : {-half,half}) {
          dcomplex M(0);
          for (uint si : {0,1,2,3}) {
            M += cutoff * ubar(laN) * vertex1hNgamma(g1,spinParity,pR,si,-k) * u(laR) * eps(si,lag) * sign_(si);
          }
          MSQR += real(M*conj(M));
        }
      }
    }
  } else if (spin==3*half) {
    u_3h u(pR);
    for (halfint laN : {-half,half}) {
      for (halfint lag : {-_1,_0,_1}) {
        for (halfint laR(-3*half); laR<=3*half; laR+=_1) {
          dcomplex M(0);
          for (uint mu : {0,1,2,3}) {
            for (uint si : {0,1,2,3}) {
              M += cutoff * ubar(laN) * vertex3hNgamma(g1,g2,g3,spinParity,mu,pR,si,-k) * u(mu,laR) * eps(si,lag) * sign_(mu)*sign_(si);
            }
          }
          MSQR += real(M*conj(M));
        }
      }
    }
  } else if (spin==5*half) {
    u_5h u(pR);
    for (halfint laN : {-half,half}) {
      for (halfint lag : {-_1,_0,_1}) {
        for (halfint laR(-5*half); laR<=5*half; laR+=_1) {
          dcomplex M(0);
          for (uint mu : {0,1,2,3}) {
            for (uint nu : {0,1,2,3}) {
              for (uint si : {0,1,2,3}) {
                M += cutoff * ubar(laN) * vertex5hNgamma(g1,g2,g3,spinParity,mu,nu,pR,si,-k) * u(mu,nu,laR) * eps(si,lag) * sign_(mu)*sign_(nu)*sign_(si);
              }
            }
          }
          MSQR += real(M*conj(M));
        }
      }
    }
  }

  if (VMD2) { MSQR *= e*e/(grho*grho); }

  return 1./(8.*pi_) * pabs/(mr*mr) * isofac/npol * MSQR;
}



double Gamma_sigma(double m) { ///< Mass dependence of total width = width of sigma -> pi + pi
  const double delta(0.3*GeV);
  const double m0 = Config::get<double>("sigma.mass");
  const double Gamma0 = Config::get<double>("sigma.width");
  const double mpi_pm = Config::get<double>("pi_pm.mass");
  const double mpi_0 = Config::get<double>("pi_0.mass");
  //  const double mpi = (2.*mpi_pm + mpi_0)/3.;
  const double mpi = mpi_pm;
  double x = m*m - 4.*mpi*mpi;
  if (x<=0) { return 0; }
  double q = sqrt(x)/2.;  // pion momentum
  double q0 = sqrt(m0*m0 - 4.*mpi*mpi)/2.; // pion momentum at the pole
  //  cerr << "q = " << q << endl
  //       << "q0 = " << q0 << endl;
  //  return Gamma0 * m0/m * pow(q/q0,3) * (q0*q0 + delta*delta)/(q*q + delta*delta);
  return Gamma0 * q/q0 * m0/m * pow((m*m/2-mpi*mpi)/(m0*m0/2-mpi*mpi),2);
}

double dGR_Nsigma_dm(string res, double mr, double mm) { // mm: mass of sigma
  const double gsi = Config::get<double>(res+".gsi");
  const halfint spin = Config::get<halfint>(res+".spin");
  const int parity = Config::get<int>(res+".parity");
  const halfint spinParity = spin*parity;
  const halfint isospin = Config::get<halfint>(res+".isospin");
  const double mn = Config::get<double>("Nucleon.mass");
  const double msi = Config::get<double>("sigma.mass");

  double npol = 2.*habs(spin) + 1.;
  double isofac = 1.;

  double pabs = sqrt(lambda(mr*mr,mn*mn,mm*mm))/(2.*mr);
  double EN = (mr*mr + mn*mn - mm*mm)/(2.*mr);
  double Esi = (mr*mr - mn*mn + mm*mm)/(2.*mr);
  FourVector pR(mr,0,0,0);
  FourVector pN(EN,0,0,-pabs);
  FourVector k(Esi,0,0,pabs);

  ubar_1h ubar(pN);

  double MSQR(0);

  if (spin==half) {
    u_1h u(pR);
    for (halfint laN : {-half,half}) {
      for (halfint laR : {-half,half}) {
	dcomplex M = ubar(laN) * vertex1hNsi(gsi,spinParity) * u(laR);
	MSQR += real(M*conj(M));
      }
    }
  } else if (spin==3*half) {
    u_3h u(pR);
    for (halfint laN : {-half,half}) {
      for (halfint laR(-3*half); laR<=3*half; laR+=_1) {
	dcomplex M(0);
	for (uint mu : {0,1,2,3}) {
	  M += ubar(laN) * vertex3hNsi(gsi,spinParity,mu,pR,-k) * u(mu,laR) * sign_(mu);
	}
	MSQR += real(M*conj(M));
      }
    }
  } else if (spin==5*half) {
    u_5h u(pR);
    for (halfint laN : {-half,half}) {
      for (halfint laR(-5*half); laR<=5*half; laR+=_1) {
	dcomplex M(0);
	for (uint mu : {0,1,2,3}) {
	  for (uint nu : {0,1,2,3}) {
	    M += ubar(laN) * vertex5hNsi(gsi,spinParity,mu,nu,pR,-k) * u(mu,nu,laR) * sign_(mu)*sign_(nu);
	  }
	}
	MSQR += real(M*conj(M));
      }
    }
  }

  double Gamma = Gamma_sigma(mm);

  //  return 1./(4.*pi_*pi_) * pabs/(mr*mr) * isofac/npol * MSQR * (mrho*mrho * Gamma)/(POW<2>(mrho*mrho-mm*mm) + mm*mm*Gamma*Gamma);
  return 1./(4.*pi_*pi_) * pabs/(mr*mr) * isofac/npol * MSQR * (mm*msi * Gamma)/(POW<2>(msi*msi-mm*mm) + msi*msi*Gamma*Gamma);
}

double GammaR_Nsigma(string res, double mr) {
  const double mn = Config::get<double>("Nucleon.mass");
  const double mpi = Config::get<double>("pi_pm.mass");
  const double mmin(2.*mpi);
  const double mmax(mr-mn);
  const int N(100);
  const double dm = (mmax-mmin)/N;
  double sum(0);
  for (double m = mmin+dm/2.; m<mmax; m+=dm) {
    sum += dGR_Nsigma_dm(res, mr, m);
  }
  return sum*dm;
}

/**
   \todo This is rather buggy. Not even the pairing of Lorentz indices is OK for spin-5/2.
 */
// double dGR_Dpi_dmdcosth(string res, double mr, double mm, double costh) { // mr: res mass, mm: running Delta mass
//   const double gDpi = Config::get<double>(res+".gDpi");
//   const halfint spin = Config::get<halfint>(res+".spin");
//   const int parity = Config::get<int>(res+".parity");
//   const halfint spinParity = spin*parity;
//   const halfint isospin = Config::get<halfint>(res+".isospin");
//   const double mn = Config::get<double>("Nucleon.mass");
//   const double mD = Config::get<double>("D1232.mass");
//   const double GD = Config::get<double>("D1232.width");
//   const double gDNpi = Config::get<double>("D1232.g0");
//   const double mpi = Config::get<double>("pi_pm.mass");

//   double npol = 2.*habs(spin) + 1.;
//   double isofac = 1.;

//   double pabs = sqrt(lambda(mr*mr,mpi*mpi,mm*mm))/(2.*mr);
//   double Epi = (mr*mr + mpi*mpi - mm*mm)/(2.*mr);
//   double ED  = (mr*mr - mpi*mpi + mm*mm)/(2.*mr);
//   FourVector pR(mr,0,0,0);
//   FourVector ppi1(Epi,0,0,-pabs);
//   FourVector pD(ED,0,0,pabs);

//   double pabs23 = sqrt(lambda(mm*mm,mpi*mpi,mn*mn))/(2.*mm);
//   double sinth = sqrt(1.-costh*costh);
//   ThreeVector pp(pabs23*sinth,0,pabs23*costh);
//   double Epi2 = (mm*mm + mpi*mpi - mn*mn)/(2.*mm);
//   double EN   = (mm*mm - mpi*mpi + mn*mn)/(2.*mm);
//   FourVector pN_tilde(EN,pp);
//   FourVector ppi2_tilde(Epi2,-pp);

//   double ch = ED/mm;
//   double sh = pabs/mm; // both indices are upper !!!
//   FourTensor Boost(ch, 0,  0, sh,
//                    0,  1,  0,  0,
//                    0,  0,  1,  0,
//                    sh, 0,  0, ch );

//   FourVector pN = Boost*pN_tilde;
//   FourVector ppi2 = Boost*ppi2_tilde;

//   // checks:
//   // FourVector pD_tilde(mm,0,0,0);
//   // FourVector pD_new = Boost*pD_tilde;
//   // PR(pD);
//   // PR(pD_new);
//   // PR(pN);
//   // PR(ppi2);
//   // PR(pN+ppi2);
//   // PR(ppi1+pN+ppi2);
  
//   const dcomplex BWD = BW(mD*mD,mD,GD); // Breit-Wigner of Delta
  
//   ubar_1h ubar(pN);

//   double MSQR(0);

//   if (spin==half) {
//     u_1h u(pR);
//     for (halfint laN : {-half,half}) {
//       for (halfint laR : {-half,half}) {
// 	for (uint mu : {0,1,2,3}) {
// 	  for (uint nu : {0,1,2,3}) {
// 	    dcomplex M = ubar(laN) * vertex3hNpi(gDNpi,3*half,mu,pD,-ppi2) * 
// 	      i_ * (gamma_(pD)+gamma_unit*mD) * BWD *
// 	      P3h(pD,mD,mu,nu) * sign_(mu) * sign_(nu) *
// 	      vertex3hNpi(gDpi,3*spinParity,nu,-pD,-ppi1) *
// 	      u(laR);
// 	    MSQR += real(M*conj(M));
// 	  }
// 	}
//       }
//     }
//   } else if (spin==3*half) {
//     u_3h u(pR);
//     for (halfint laN : {-half,half}) {
//       for (halfint laR(-3*half); laR<=3*half; laR+=_1) {
// 	dcomplex M(0);
// 	for (uint mu : {0,1,2,3}) {
// 	  for (uint nu : {0,1,2,3}) {
// 	    for (uint ro : {0,1,2,3}) {
// 	      M += ubar(laN) * vertex3hNpi(gDNpi,3*half,mu,pD,-ppi2) * 
// 		i_ * (gamma_(pD)+gamma_unit*mD) * BWD *
// 		P3h(pD,mD,mu,nu) * sign_(mu) * sign_(nu) *
// 		vertex3h3hpi(gDpi,3*half,nu,-pD,spinParity,ro,pR,-ppi1) *
// 		u(mu,laR);
// 	    }
// 	  }
// 	}
// 	MSQR += real(M*conj(M));
//       }
//     }
//   // } else if (spin==5*half or spin==-5*half) {
//   //   u_5h u(pR);
//   //   for (halfint laN : {-half,half}) {
//   //     for (halfint laR(-5*half); laR<=5*half; laR+=1) {
//   // 	dcomplex M(0);
//   // 	for (uint mu : {0,1,2,3}) {
//   // 	  for (uint nu : {0,1,2,3}) {
//   // 	    M += ubar(laN) * vertex5hNsi(gsi,spin,mu,nu,pR,-k) * u(mu,nu,laR);
//   // 	    //M += ubar(laN) * vertex5hNsi(gsi,spin,mu,nu,pR,-k) * u(mu,nu,laR) * sign_(mu)*sign_(nu);
//   // 	  }
//   // 	}
//   // 	MSQR += real(M*conj(M));
//   //     }
//   //   }
//   }

//   return 1./(16.*POW<3>(2.*pi_)) * pabs*pabs23/(mr*mr*mm) * isofac/npol * MSQR;
// }


// double dGR_Dpi_dm(string res, double mr, double mm) {
//   const int N(50);
//   double sum(0);
//   double dcosth = 2./N;
//   for (double costh(-1. + dcosth/2.); costh<1.; costh+=dcosth) {
//     sum += dGR_Dpi_dmdcosth(res,mr,mm,costh);
//   }
//   return sum * dcosth;
// }

// double GammaR_Dpi(string res, double mr) {
//   const double mn = Config::get<double>("Nucleon.mass");
//   const double mpi = Config::get<double>("pi_pm.mass");
//   const double mmin(mn + mpi);
//   const double mmax(mr - mpi);
//   const int N(50);
//   const double dm = (mmax-mmin)/N;
//   double sum(0);
//   for (double m = mmin+dm/2.; m<mmax; m+=dm) {
//     sum += dGR_Dpi_dm(res,mr,m);
//   }
//   return sum*dm;
// }

double GammaR_tot(std::string res, double mr) {
  return GammaR_Npi(res,mr)
    + GammaR_Nrho(res,mr)
    + GammaR_Nsigma(res,mr)
    //+ GammaR_Dpi(res,mr)
    ;
}

bindResonance::bindResonance(std::function<double(std::string,double)> func, std::string resonance) :
  func(func),
  resonance(resonance),
  scaleFactor(getScaleFactor(func,resonance)) {
}

double bindResonance::operator()(double x) const {
  return scaleFactor*func(resonance,x);
}

double bindResonance::getScaleFactor(std::function<double(std::string,double)> func, std::string resonance) {
  double scaleFactor(1.);
  if (Config::exists("rescale")) {
    double rmass = Config::get<double>(resonance+".mass");
    double rwidth = Config::get<double>(resonance+".width");
    double act_width = func(resonance,rmass);
    scaleFactor = rwidth/act_width;
  }
  return scaleFactor;
}
