#include "resdec.hpp"

#include <iostream>
#include <iomanip>

#include "Config.hpp"

using namespace std;

int main(int argc, char** argv) {
  Config::load(argc,argv);
  /*
  dGR_Nrho_dm("N1440", 1.44, 0.5);
  exit(0);
  //*/
  /*
  dGR_Nrho_dm("N1520", 1.52, 0.5041);
  exit(0);
  //*/
  /*
  dGR_Nrho_dm("N1680", 1.68, 0.5);
  exit(0);
  //*/
  double msi = Config::get<double>("sigma.mass");
  double Gamma_sipipi = Config::get<double>("sigma.width");
  double gsipipi = Config::get<double>("gsipipi");
  double tmp_sipipi = Gsigma_pipi(gsipipi,msi);
  double gsipipi_phys = gsipipi*sqrt(Gamma_sipipi/tmp_sipipi);
  cout << "sigma" << endl
       << "      gsipipi:  " << setw(10) << gsipipi << " -> " << setw(10) << gsipipi_phys
       << "    ( " << setw(10) << tmp_sipipi << " -> " << setw(10) << Gamma_sipipi << " )" << endl;

  for (string res : {"D1232","N1440","N1520","N1535","N1650","N1675","N1680","N1700","N1710","D1600","D1620"}) {
    cout << "res: " << res << endl;
    const double mr = Config::get<double>(res+".mass");
    double tmp_ngamma = GammaR_Ngamma(res,mr);
    double tmp_ngamma_VMD2 = GammaR_Ngamma(res,mr,"VMD2");
    double tmp_pi = GammaR_Npi(res,mr);
    double tmp_rho = GammaR_Nrho(res,mr);
    double tmp_si = GammaR_Nsigma(res,mr);
    //double tmp_Dpi = GammaR_Dpi(res,mr);
    double width = Config::get<double>(res+".width");
    double BR_ngamma = Config::get<double>(res+".Bngamma");
    double BR_pi = Config::get<double>(res+".BNpi");
    double BR_rho = Config::get<double>(res+".BNrho");
    double BR_si = Config::get<double>(res+".BNsi");
    //double BR_Dpi = Config::get<double>(res+".BDpi");
    double Gamma_ngamma = width*BR_ngamma;
    double Gamma_pi = width*BR_pi;
    double Gamma_ro = width*BR_rho;
    double Gamma_si = width*BR_si;
    //double Gamma_Dpi = width*BR_Dpi;
    double gngamma = Config::get<double>(res+".gngamma");
    double gngamma_phys = gngamma*sqrt(Gamma_ngamma/tmp_ngamma);
    double gVMD2 = Config::get<double>(res+".gVMD2");
    double gVMD2_phys = gVMD2*sqrt(Gamma_ngamma/tmp_ngamma_VMD2);
    double g0 = Config::get<double>(res+".g0");
    double g0_phys = g0*sqrt(Gamma_pi/tmp_pi);
    double g1 = Config::get<double>(res+".g1");
    double g1_phys = g1*sqrt(Gamma_ro/tmp_rho);
    double gsi = Config::get<double>(res+".gsi");
    double gsi_phys = gsi*sqrt(Gamma_si/tmp_si);
    //double gDpi = Config::get<double>(res+".gDpi");
    //double gDpi_phys = gDpi*sqrt(Gamma_Dpi/tmp_Dpi);
    
    cout << res << endl
	 << "  gngamma:  " << setw(10) << gngamma << " -> " << setw(10) << gngamma_phys
         << "    ( " << setw(10) << tmp_ngamma << " -> " << setw(10) << Gamma_ngamma << " )" << endl
	 << "    gVMD2:  " << setw(10) << gVMD2 << " -> " << setw(10) << gVMD2_phys
         << "    ( " << setw(10) << tmp_ngamma_VMD2 << " -> " << setw(10) << Gamma_ngamma << " )" << endl
	 << "       g0:  " << setw(10) << g0 << " -> " << setw(10) << g0_phys
         << "    ( " << setw(10) << tmp_pi << " -> " << setw(10) << Gamma_pi << " )" << endl
         << "       g1:  " << setw(10) << g1 << " -> " << setw(10) << g1_phys
         << "    ( " << setw(10) << tmp_rho << " -> " << setw(10) << Gamma_ro << " )" << endl
         << "      gsi:  " << setw(10) << gsi << " -> " << setw(10) << gsi_phys
         << "    ( " << setw(10) << tmp_si << " -> " << setw(10) << Gamma_si << " )" << endl;
         // << "     gDpi:  " << setw(10) << gDpi << " -> " << setw(10) << gDpi_phys
         // << "    ( " << setw(10) << tmp_Dpi << " -> " << setw(10) << Gamma_Dpi << " )" << endl;
  }

  //  cout << "rho and sigma width parametrization should be:   Gamma0 * pow(q/q0,3) * m0/m    !!!!!" << endl;
  //  cout << "isospin factors fro Delta + pion decay are missing!" << endl;
  return 1;
}
