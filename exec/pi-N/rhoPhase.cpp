#include <iostream>
#include <iomanip>
#include "utils.hpp"
#include "Gamma.hpp"
#include "Config.hpp"

using namespace std;

int main(int argc, char** argv) {

    Config::load(argc,argv);

    double Mmin(0);
    if (Config::exists("Mmin")) Mmin = Config::get<double>("Mmin");
    double Mmax(1.);
    if (Config::exists("Mmax")) Mmin = Config::get<double>("Mmax");
    double dM(0.01);
    if (Config::exists("dM")) dM = Config::get<double>("dM");
    double mrho = Config::get<double>("rho.mass");

    for (double M(Mmin); M<Mmax; M+=dM) {
        dcomplex BWrho = -BW(M*M,mrho,Gamma_rho(M));
        //dcomplex BWrho = -1./(M*M - mrho*mrho + i_*M*Gamma_rho(M));
        double phi = arg(BWrho);
        cout << setw(12) << M << setw(14) << phi*180./pi_ << endl;
    }
    return 0;
}