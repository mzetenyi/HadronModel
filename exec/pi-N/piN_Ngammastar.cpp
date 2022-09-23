#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "piN_Ngammastar.hpp"
#include "Gamma.hpp"
#include "BornTerms.hpp"
#include "MultiArray.hpp"
#include "wavefunc.hpp"
#include "RandomNumberGenerator.hpp"
#include "Config.hpp"
#include "Histogram.hpp"
#include "cutoffs.hpp"

using namespace std;

vector<string> piN_Ngammastar::determineResonances() const {
  //  cerr << "in vector<string> piN_Ngammastar::determineResonances() const" << endl;
  vector<string> resonances;
  if (Config::exists("D1232")) resonances.push_back("D1232");
  if (Config::exists("N1440")) resonances.push_back("N1440");
  if (Config::exists("N1535")) resonances.push_back("N1535");
  if (Config::exists("N1520")) resonances.push_back("N1520");
  if (Config::exists("N1650")) resonances.push_back("N1650");
  if (Config::exists("N1675")) resonances.push_back("N1675");
  if (Config::exists("N1680")) resonances.push_back("N1680");
  if (Config::exists("N1710")) resonances.push_back("N1710");
  if (Config::exists("D1600")) resonances.push_back("D1600");
  if (Config::exists("D1620")) resonances.push_back("D1620");
  //  cerr << "return from vector<string> piN_Ngammastar::determineResonances() const" << endl;
  return resonances;
}

vector<string> piN_Ngammastar::determineChannels() const {
  //  cerr << "in vector<string> piN_Ngammastar::determineChannels() const" << endl;
  vector<string> channels;
  if (Config::exists("Born")) {
    channels.push_back("Born");
  } else {
    if (Config::exists("Born_s")) channels.push_back("Born_s");
    if (Config::exists("Born_u")) channels.push_back("Born_u");
    if (Config::exists("Born_t")) channels.push_back("Born_t");
    if (Config::exists("Born_c")) channels.push_back("Born_c");
    if (Config::exists("Born_corr")) channels.push_back("Born_corr");
  }
  for (string res : resonances) {
    channels.push_back(res);
  }
  //  cerr << "return from vector<string> piN_Ngammastar::determineChannels() const" << endl;
  return channels;
}

vector<string> piN_Ngammastar::determineEMchannels() const {
  //  cerr << "in vector<string> piN_Ngammastar::determineEMchannels() const" << endl;
  vector<string> EMchannels;
  bool gamma = Config::exists("gamma");
  bool rho = (Config::exists("rho") or Config::exists("rho2"));
  if (not (gamma or rho)) gamma = rho = true;
  bool Bgamma = Config::exists("Bgamma");
  bool Brho = (Config::exists("Brho") or Config::exists("Brho2"));
  if (gamma or Bgamma) EMchannels.push_back("gamma");
  if (rho or Brho) EMchannels.push_back("rho");
  //  cerr << "return from vector<string> piN_Ngammastar::determineEMchannels() const" << endl;
  return EMchannels;
}

// vector<NamedIndex> getChannelIndices() const {
//   vector<NamedIndex> ChInd;
//   ChInd.push_back(NamedIndex(channels));
//   ChInd.push_back(NamedIndex(EMchannels));
//   return ChInd;
// }

// vector<string> piN_Ngammastar::getCouplings() const {
//   vector<string> couplings;
//   for (string res : resonances) {
//     for (string decCh : {"Npi","Ngamma","Nrho"}) {
//       couplings.push_back(res+decCh);
//     }
//   }
//   return couplings;
// }

MultiArray<double> piN_Ngammastar::getRelative_errors() const {
  MultiArray<double> relative_errors({NIch,NIdec});
  for (string res : channels) {
    double u_Npi = 0;
    double u_Ngamma = 0;
    double u_Nrho = 0;
    if (Config::exists(res+".u_g0")) u_Npi = Config::get<double>(res+".u_g0");
    if (Config::exists(res+".u_gngamma")) u_Ngamma = Config::get<double>(res+".u_gngamma");
    if (Config::exists(res+".u_g1")) u_Nrho = Config::get<double>(res+".u_g1");
    relative_errors(NIch(res),NIdec("pi")) = u_Npi;
    relative_errors(NIch(res),NIdec("gamma")) = u_Ngamma;
    relative_errors(NIch(res),NIdec("rho")) = u_Nrho;
    if (Config::exists("rho2") and Config::exists("VMD2")) relative_errors(NIch(res),NIdec("rho")) = u_Ngamma;
  }
  return relative_errors;
}

piN_Ngammastar::piN_Ngammastar(double srt, double M) : 
  srt(srt),
  M(M),
  resonances(determineResonances()),
  channels(determineChannels()),
  EMchannels(determineEMchannels()),
  DECchannels({"pi","gamma","rho"}),
  NIch(channels),
  NIem(EMchannels),
  NIdec(DECchannels),
  relative_errors(getRelative_errors()),
  Bgamma(Config::exists("Bgamma")),
  Brho(Config::exists("Brho")),
  Brho2(Config::exists("Brho2")),
  gamma(Config::exists("gamma")),
  rho(Config::exists("rho")),
  rho2(Config::exists("rho2")),
  VMD2(Config::exists("VMD2")),
  vRho(vectorRho()),
  tRho(tensorRho()),
  vGamma(vectorGamma()),
  tGamma(tensorGamma()),
  vOmega(vectorOmega()),
  tOmega(tensorOmega()),
  sch(Config::exists("sch")),
  uch(Config::exists("uch")),
  Born(Config::exists("Born")),
  Born_s(Config::exists("Born_s")),
  Born_u(Config::exists("Born_u")),
  Born_t(Config::exists("Born_t")),
  Born_c(Config::exists("Born_c")),
  Born_corr(Config::exists("Born_corr")) {
  //  cerr << "in piN_Ngammastar constructor" << endl;
  if (not(gamma or rho or rho2)) gamma = rho = true;
  if (rho2) gamma = rho = false;
  if (Brho2) Bgamma = Brho = false;
  if (not(Bgamma or Brho or Brho2)) {
    Bgamma = gamma;
    Brho = rho;
    Brho2 = rho2;
  }
  if (not(sch or uch)) sch = uch = true;
  s = srt*srt;
  M2 = M*M;
  mpi_pm = Config::get<double>("pi_pm.mass");
  mpi_pm2 = mpi_pm*mpi_pm;
  mpi_0 = Config::get<double>("pi_0.mass");
  mpi = mpi_pm;
  //mpi = (2.*mpi_pm + mpi_0)/3.;
  mpi2 = mpi*mpi;
  mrho = Config::get<double>("rho.mass");
  mn = Config::get<double>("Nucleon.mass");
  mn2 = mn*mn;
  if (srt<mn+mpi) { cerr << "sqrt(s) is below threshold of pi+N" << endl; exit(0); }
  if (srt<mn+M) { cerr << "sqrt(s) is below threshold of N + dilep of mass M" << endl; exit(0); }
  EN_in = (s-mpi2+mn2)/(2.*srt);
  Epi_in = (s+mpi2-mn2)/(2.*srt);
  pin_abs = sqrt(lambda(s,mpi2,mn2))/(2.*srt);
  //  cerr << "incoming threemomentum: " << pin_abs << endl;
  p1 = FourVector(EN_in,0,0,-pin_abs);
  q = FourVector(Epi_in,0,0,pin_abs);
  //P = p1+p2;
  EN_out = (s-M2+mn2)/(2.*srt);
  k0 = (s+M2-mn2)/(2.*srt);
  pout_abs = sqrt(lambda(s,M2,mn2))/(2.*srt);
  //  cerr << "outgoing threemomentum: " << pout_abs << endl;
  const double grho_tilde = Config::get<double>("grho_tilde");
  const double grho = Config::get<double>("grho");
  F_rho = grho_tilde/grho * (M2)/(M2 - mrho*mrho + i_*sqrt(M2)*Gamma_rho(M));
  if (Config::exists("noFrho")) {
    cerr << "F_rho is set to 1 via the option \"noFrho\" !!!!   Orig. value would be: " << F_rho << endl;
    //cerr << "|F_rho|^2 = " << F_rho*conj(F_rho)*(e*e)/(grho_tilde*grho_tilde) << endl;
    F_rho = 1.;
  }
  //cerr << "  e = " << e << endl;
  // cerr << "  grho = " << grho << endl;
  // cerr << "  grho_tilde = " << grho_tilde << endl;
  //F_rho = -e/grho * (M2)/(M2 - mrho*mrho + i_*mrho*Gamma_rho(M));
  //cerr << "F_rho = " << F_rho << endl;
  // cerr << "F_rho^2 = " << F_rho*conj(F_rho) << endl;
  //  cerr << "end of constructor" << endl;
}

double piN_Ngammastar::u_Mandelstam(double costh) const {
  double sinth = sqrt(1-costh*costh);
  ThreeVector kk = pout_abs*ThreeVector(sinth,0,costh);
  FourVector k = FourVector(k0,kk);
  return (p1-k)*(p1-k);
}
//*


double piN_Ngammastar::cutoff_sch(const std::string& res) const {
  const double mr = Config::get<double>(res+".mass");
  const double Gr = Config::get<double>(res+".width");
  const double BNpi = Config::get<double>(res+".BNpi");
  const int l = Config::get<int>(res+".l");
  const halfint spin = Config::get<halfint>(res+".spin");
  
  if (Config::exists("cutoff")) {
    if (Config::get<string>("cutoff")=="Rapp") {
      double qpi2(pin_abs*pin_abs);
      //      double qrho2(pout_abs*pout_abs);
      //      return cutoffRapp(spin,qpi2) * cutoffRapp(spin,qrho2);
      return cutoffRapp(spin,qpi2);
    } else {
      double mr2(mr*mr);
      double delta2;
      double Gamma_r = Gr*BNpi;
      if (res == "D1232") {
        delta2 = POW<2>(0.3*GeV);
      } else if (res == "N1535") {
        delta2 = POW<2>(0.5*GeV);
      } else {
        delta2 = POW<2>(mr-mn-mpi) + POW<2>(Gamma_r)/4.;
      }
      double q2(pin_abs*pin_abs);
      double q02 = lambda(mr2,mn2,mpi2)/(4.*mr2);
      return pow((q02+delta2)/(q2+delta2), (l+1.)/2.) * sqrt(mr/srt);
    }
  }
  return 1.;
}



/**
   Calculate the helicity amplitudes
*/
HelicityAmplitudes piN_Ngammastar::helicityAmplitudes(double costh) const {
  // isospin factors:
  const double isoNs = -sqrt(2.);
  const double isoNu =  sqrt(2.);
  const double isoDs = sqrt(2.)/3.;
  const double isoDu = -sqrt(2.)/3.;

  const double isoNds = sqrt(2.);
  const double isoNdu = sqrt(2.);
  const double isoDds = 1./sqrt(3.);
  const double isoDdu = -1./sqrt(3.);

  const double iso = sqrt(2.); // for Born terms

  double sinth = sqrt(1-costh*costh);
  ThreeVector kk = pout_abs*ThreeVector(sinth,0,costh);
  FourVector p2 = FourVector(EN_out,-kk);
  FourVector k = FourVector(k0,kk);

  double ch = k0/M;
  double sh = -pout_abs/M;
  FourTensor Boost = FourTensor::createUpDown( ch,  0,  0, -sh,
                                               0,   1,  0,  0,
                                               0,   0,  1,  0,
                                               -sh, 0,  0,  ch );

  FourTensor Rot = FourTensor::createUpDown( 1,    0,   0,   0,
                                             0,  costh, 0, sinth,
                                             0,    0,   1,   0,
                                             0, -sinth, 0, costh );

  FourVector ktilde(M,0,0,0);
  FourVector kprime = Boost*ktilde;
  FourVector knew = Rot*kprime;

  MultiArray<FourVector> epsRe(idx_s1);
  MultiArray<FourVector> epsIm(idx_s1);

  // photon wave vectors:
  FourVector epsRe0(0,0,0,1.);
  FourVector epsIm0(0,0,0,0);
  epsRe(_0) = Rot*Boost*epsRe0;
  epsIm(_0) = Rot*Boost*epsIm0;

  FourVector epsRep(0,-1./sqrt(2.),0,0);
  FourVector epsImp(0,0,-1./sqrt(2.),0);
  epsRe(_1) = Rot*Boost*epsRep;
  epsIm(_1) = Rot*Boost*epsImp;

  FourVector epsRem(0,1./sqrt(2.),0,0);
  FourVector epsImm(0,0,-1./sqrt(2.),0);
  epsRe(-_1) = Rot*Boost*epsRem;
  epsIm(-_1) = Rot*Boost*epsImm;

  double uu = (p1-k)*(p1-k);
  double tt = (p1-p2)*(p1-p2);

  // form factors for Born terms:
  double F1 = 1.;
  double F2 = 1.;
  double F3 = 1.;
  double Ftilde = 1.;

  if (not Config::exists("noBornCutoff")) {
    double LBorn = 0.63; // from Zetenyi, Wolf PRC 2012
    if (Config::exists("LBorn")) {
      LBorn = Config::get<double>("LBorn");
    }

    F1 = 1./(1 + POW<2>(s-mn2)/POW<4>(LBorn));
    F2 = 1./(1 + POW<2>(uu-mn2)/POW<4>(LBorn));
    F3 = 1./(1 + POW<2>(tt-mpi2)/POW<4>(LBorn));
    Ftilde = F1 + F2 + F3 - F1*F2 - F1*F3 - F2*F3 + F1*F2*F3;
  }

  MultiArray<dcomplex> Mhad(NIch.getIdx(),NIem.getIdx(),idx_s1h,idx_s1h,idx_lor);

  ubar_ ubar(half,p2);
  u_ u(half,p1);

  MultiArray<DiracMatrix> M_(idx(_1,_5),idx_lor);
  for (halfint mu(0); mu<4; mu++) {
    M_(_1,mu) = gamma5_*(gamma_(mu)*gamma_(k) - k(mu)*gamma_unit);
    M_(_2,mu) = gamma5_/2.*((2.*p1(mu)-k(mu))*(2.*(q*k)-M2) - (2.*q(mu)-k(mu))*(2.*(p1*k)-M2));
    M_(_3,mu) = gamma5_/2.*(gamma_(mu)*(2.*(p2*k)+M2) - (2.*p2(mu)+k(mu))*gamma_(k));
    M_(_4,mu) = gamma5_/2.*(gamma_(mu)*(2.*(p1*k)-M2) - (2.*p1(mu)-k(mu))*gamma_(k));
    M_(_5,mu) = gamma5_/2.*((2.*p2(mu)+k(mu))*(2.*(q*k)-M2) - (2.*q(mu)-k(mu))*(2.*(p2*k)+M2));
  }

  dcomplex Agamma[6];
  dcomplex Arho[6];
  const double kan = Config::get<double>("ka_nngamma");
  const double kap = Config::get<double>("ka_ppgamma");
  const double kar = Config::get<double>("ka_NNrho");
  const double f_NNpi = Config::get<double>("f_NNpi");
  const double grho_tilde = Config::get<double>("grho_tilde");

  for (uint i(1); i<=5; i++) {
    Agamma[i] = Arho[i] = 0.;
  }

  if (Bgamma) {
    if (vGamma) {
      Agamma[1] += 2.*mn*F2/(uu-mn2);
      Agamma[2] += - 4.*mn*Ftilde/(uu-mn2)/(tt-mpi2);
    }

    if (tGamma) {
      if (Config::exists("onlyS")) {
        Agamma[1] += - 2.*mn*kan*F1/(s-mn2) - kan*F1/(2.*mn);
        Agamma[3] += 2.*kan*F1/(s-mn2);
      } else if (Config::exists("onlyU")) {
        Agamma[1] += - 2.*mn*kap*F2/(uu-mn2) - kap*F2/(2.*mn);
        Agamma[4] += 2.*kap*F2/(uu-mn2);
      } else {
        Agamma[1] += - 2.*mn*(kan*F1/(s-mn2) + kap*F2/(uu-mn2)) - (kan*F1 + kap*F2)/(2.*mn);
        Agamma[3] += 2.*kan*F1/(s-mn2);
        Agamma[4] += 2.*kap*F2/(uu-mn2);
      }
    }
  }

  if (Brho) {
    if (vRho) {
      Arho[1] += mn*(F2/(uu-mn2)-F1/(s-mn2)) * F_rho;
      Arho[2] += - 2.*mn*Ftilde/(uu-mn2)/(tt-mpi2) * F_rho;
      Arho[5] += 2.*mn*Ftilde/(s-mn2)/(tt-mpi2) * F_rho;
    }

    if (tRho) {
      Arho[1] += (- 2.*mn*(-kar*F1/(s-mn2) + kar*F2/(uu-mn2)) - (-kar*F1 + kar*F2)/(2.*mn)) * F_rho/2.;
      Arho[3] += -2.*kar*F1/(s-mn2) * F_rho/2.;
      Arho[4] += 2.*kar*F2/(uu-mn2) * F_rho/2.;
    }
  }

  if (Brho2) {
    if (vRho) {
      Arho[1] += mn*(F2/(uu-mn2)-F1/(s-mn2)) * (mrho*mrho/M2)*F_rho;
      Arho[2] += - 2.*mn*Ftilde/(uu-mn2)/(tt-mpi2) * (mrho*mrho/M2)*F_rho;
      Arho[5] += 2.*mn*Ftilde/(s-mn2)/(tt-mpi2) * (mrho*mrho/M2)*F_rho;
    }

    if (tRho) {
      Arho[1] += (- 2.*mn*(-kar*F1/(s-mn2) + kar*F2/(uu-mn2)) - (-kar*F1 + kar*F2)/(2.*mn)) * (mrho*mrho/M2)*F_rho/2.;
      Arho[3] += -2.*kar*F1/(s-mn2) * (mrho*mrho/M2)*F_rho/2.;
      Arho[4] += 2.*kar*F2/(uu-mn2) * (mrho*mrho/M2)*F_rho/2.;
    }
  }
  
  dcomplex Born_rhofac(1.);
  if (Config::exists("Bornrhophase")) {
    double phase = Config::get<double>("Bornrhophase"); // phase in degrees
    Born_rhofac = cos(phase*pi_/180.) + i_*sin(phase*pi_/180.);
  }
  
  for (int j(1); j<=5; j++) {
    Agamma[j] *= e*f_NNpi*sqrt(2)/mpi;
    Arho[j]   *= e*f_NNpi*sqrt(2)/mpi * Born_rhofac;
  }

  for (halfint la1 : {-half,half}) {
    for (halfint la2 : {-half,half}) {
      for (halfint mu(0); mu<4; mu++) {
        for (string ch : channels) {
          for (string emch : EMchannels) {
            Mhad(NIch(ch),NIem(emch),la1,la2,mu) = 0;
          }
        }
        //------------  Born terms  ------------------------------
        if (Born) {
          for (int j(1); j<=5; j++) {
            if (Bgamma)        Mhad(NIch("Born"),NIem("gamma"),la1,la2,mu) += ubar(0,la2) * Agamma[j] * M_(_1*j,mu) * u(0,la1);
            if (Brho or Brho2) Mhad(NIch("Born"),NIem("rho"),la1,la2,mu)   += ubar(0,la2) *   Arho[j] * M_(_1*j,mu) * u(0,la1);
          }
        }
        //------------  Born terms separately --------------------
        if (Born_s) { // s-channel
          if (Bgamma) {
            Mhad(NIch("Born_s"),NIem("gamma"),la1,la2,mu) += F1 * ubar(0,la2) * vertexNNgamma(mu, 0, -k) *
              fprop(p2+k,mn) * vertexNNpi(q,1,-1) * u(0,la1);
          }
          if (Brho) {
            Mhad(NIch("Born_s"),NIem("rho"),la1,la2,mu) += F1 * ubar(0,la2) * vertexNNrho(mu, -k,0,0) *
              fprop(p2+k,mn) * vertexNNpi(q,1,-1) * u(0,la1) * (-e/grho_tilde)*F_rho;
          }
        }
        if (Born_u) { // u-channel
          if (Bgamma) {
            Mhad(NIch("Born_u"),NIem("gamma"),la1,la2,mu) += F2 * ubar(0,la2) * vertexNNpi(q,1,-1) *
              fprop(p1-k,mn) * vertexNNgamma(mu, 1, -k) * u(0,la1);
          }
          if (Brho) {
            Mhad(NIch("Born_u"),NIem("rho"),la1,la2,mu) += F2 * ubar(0,la2) * vertexNNpi(q,1,-1) *
              fprop(p1-k,mn) * vertexNNrho(mu, -k,1,0) * u(0,la1) * (-e/grho_tilde)*F_rho;
          }
        }
        if (Born_t) { // t-channel
          if (Bgamma) {
            Mhad(NIch("Born_t"),NIem("gamma"),la1,la2,mu) += F3 * ubar(0,la2) * vertexNNpi(-(k-q),1,-1) * u(0,la1) *
              sprop(k-q,mpi) * vertexgammapipi(mu, k-q, q);
          }
          if (Brho) {
            Mhad(NIch("Born_t"),NIem("rho"),la1,la2,mu) += F3 * ubar(0,la2) * vertexNNpi(-(k-q),1,-1) * u(0,la1) *
              sprop(k-q,mpi) * vertexrhopipi(mu, k-q, q) * (-e/grho_tilde)*F_rho;
          }
        }
        if (Born_c) { // contact term
          if (Bgamma) {
            Mhad(NIch("Born_c"),NIem("gamma"),la1,la2,mu) += iso * ubar(0,la2) * vertexNNgammapi(mu) * u(0,la1);
          }
          if (Brho) {
            Mhad(NIch("Born_c"),NIem("rho"),la1,la2,mu) += ubar(0,la2) * vertexNNrhopi(mu,1,0,-1) * u(0,la1) * (-e/grho_tilde)*F_rho;
          }
        }
        if (Born_corr) { // correction term for gauge invariance
          if (Bgamma) {
            Mhad(NIch("Born_corr"),NIem("gamma"),la1,la2,mu) += -e*f_NNpi*iso/mpi * ubar(0,la2) * 
              (
               - 2.*mn*gamma5_ * (
                                  (F1-Ftilde)*(k(mu)-2.*p1(mu))/(s-mn2)
                                  - (F3-Ftilde)*(k(mu)-2.*q(mu))/(tt-mpi2)
                                  )
               - gamma5_*gamma_(mu) * ((F1-Ftilde) - (F3-Ftilde))
               ) *
              u(0,la1);
          }
          if (Brho) {
            Mhad(NIch("Born_corr"),NIem("rho"),la1,la2,mu) += grho_tilde*f_NNpi/(sqrt(2.)*mpi) * 2.*mn * ubar(0,la2) * gamma5_ * 
              ((Ftilde-F2)*(2.*p1(mu)-k(mu))/(uu-mn2)
               - (Ftilde-F1)*(2.*p2(mu)+k(mu))/(s-mn2)
               - 2.*(Ftilde-F3)*(2.*q(mu)-k(mu))/(tt-mpi2)) *
              u(0,la1);
          }
        }

        // ---------- Resonance contributions ----------------------
        for (const string& res : resonances) {
          if (Config::exists(res)) {
            const halfint spin = Config::get<halfint>(res+".spin");
            const int parity = Config::get<int>(res+".parity");
            const halfint spinparity = spin*parity;
            const halfint isospin = Config::get<halfint>(res+".isospin");
            const double mr = Config::get<double>(res+".mass");
            const double Gr = Config::get<double>(res+".width");
            const double gA = Config::get<double>(res+".gngamma");
            const double g0 = Config::get<double>(res+".g0");
            // VMD2 below: use rho coupling fixed at photon point 
            const double g1 = (VMD2 ? Config::get<double>(res+".gVMD2") : Config::get<double>(res+".g1"));
            const double g2 = Config::get<double>(res+".g2");
            const double g3 = Config::get<double>(res+".g3");
            const int l = Config::get<int>(res+".l");
            const dcomplex BWs = BW(s,mr,Gamma_R(res,s,mr,Gr,l));
            const dcomplex BWu = BW(uu,mr,0);
            double iso_s(0);
            double iso_u(0);
            double iso_ds(0);
            double iso_du(0);
            if (isospin==half) {
              iso_s = isoNs;
              iso_u = isoNu;
              iso_ds = isoNds;
              iso_du = isoNdu;
            } else if (isospin==3*half) {
              iso_s = isoDs;
              iso_u = isoDu;
              iso_ds = isoDds;
              iso_du = isoDdu;
            } else {
              cerr << "unknown isospin for resonance " << res << ": " << isospin << endl;
              exit(0);
            }
            dcomplex fac = 1.;
            if (Config::exists(res+".phase")) {
              double phase = Config::get<double>(res+".phase"); // phase in degrees
              fac = cos(phase*pi_/180.) + i_*sin(phase*pi_/180.);
            }
            dcomplex rhofac(1.);
            if (Config::exists(res+".rhophase")) {
              double phase = Config::get<double>(res+".rhophase"); // phase in degrees
              rhofac = cos(phase*pi_/180.) + i_*sin(phase*pi_/180.);
            } else if (Config::exists("res_rhophase")) {
              double phase = Config::get<double>("res_rhophase"); // phase in degrees
              rhofac = cos(phase*pi_/180.) + i_*sin(phase*pi_/180.);
            }
            // cutoff:
            double cutoff = cutoff_sch(res);
            //------------  spin-1/2  ------------------------------
            if (spin == half) {
              if (gamma) {
                Mhad(NIch(res),NIem("gamma"),la1,la2,mu) += 
                  fac *
                  ubar(0,la2) *
                  ( iso_ds *
                    ( (not sch) ? gamma_null
                      : ( cutoff * vertex1hNgamma(gA,spinparity,p1+q,mu,-k) *
                          i_ * (gamma_(p1+q)+gamma_unit*mr) * BWs *
                          vertex1hNpi(g0,spinparity,q) ) )
                    + iso_du *
                    ( (not uch) ? gamma_null
                      : ( vertex1hNpi(g0,spinparity,q) *
                          i_ * (gamma_(p1-k)+gamma_unit*mr) * BWu *
                          vertex1hNgamma(gA,spinparity,-p1+k,mu,-k) )
                      )
                    ) * u(0,la1);
              }
              if (rho) {
                Mhad(NIch(res),NIem("rho"),la1,la2,mu) += 
                  fac*rhofac * (-e/grho_tilde)*F_rho *
                  ubar(0,la2) *
                  ( iso_s *
                    ( (not sch) ? gamma_null
                      : ( cutoff * vertex1hNrho(g1,spinparity,p1+q,mu,-k) *
                          i_ * (gamma_(p1+q)+gamma_unit*mr) * BWs *
                          vertex1hNpi(g0,spinparity,q) ) )
                    + iso_u *
                    ( (not uch) ? gamma_null
                      : ( vertex1hNpi(g0,spinparity,q) *
                          i_ * (gamma_(p1-k)+gamma_unit*mr) * BWu *
                          vertex1hNrho(g1,spinparity,-p1+k,mu,-k) )
                      )
                    ) * u(0,la1);
              }
              if (rho2) {
                Mhad(NIch(res),NIem("rho"),la1,la2,mu) += 
                  fac*rhofac * (-e/grho_tilde)*(mrho*mrho/M2)*F_rho * 
                  ubar(0,la2) *
                  ( iso_s *
                    ( (not sch) ? gamma_null
                      : ( cutoff * vertex1hNrho(g1,spinparity,p1+q,mu,-k) *
                          i_ * (gamma_(p1+q)+gamma_unit*mr) * BWs *
                          vertex1hNpi(g0,spinparity,q) ) )
                    + iso_u *
                    ( (not uch) ? gamma_null
                      : ( vertex1hNpi(g0,spinparity,q) *
                          i_ * (gamma_(p1-k)+gamma_unit*mr) * BWu *
                          vertex1hNrho(g1,spinparity,-p1+k,mu,-k) )
                      )
                    ) * u(0,la1);
              }
            }
            //------------  spin-3/2  ------------------------------
            if (spin == 3*half) {
              for (uint al(0); al<4; al++) {
                for (uint be(0); be<4; be++) {
                  if (gamma) {
                    Mhad(NIch(res),NIem("gamma"),la1,la2,mu) +=
                      fac * 
                      ubar(0,la2) *
                      ( iso_ds *
                        ( (not sch) ? gamma_null
                          : ( cutoff * vertex3hNgamma(gA,0,0,spinparity,al,p1+q,mu,-k) *
                              i_ * (gamma_(p1+q)+gamma_unit*mr) * BWs *
                              P3h(p1+q,mr,al,be)*sign_(al)*sign_(be) *
                              vertex3hNpi(g0,spinparity,be,-p1-q,q) ) )
                        + iso_du *
                        ( (not uch) ? gamma_null
                          : ( vertex3hNpi(g0,spinparity,be,p1-k,q) *
                              i_ * (gamma_(p1-k)+gamma_unit*mr) * BWu *
                              P3h(p1-k,mr,be,al)*sign_(al)*sign_(be) *
                              vertex3hNgamma(gA,0,0,spinparity,al,-p1+k,mu,-k) ) )
                        )
                      * u(0,la1);
                  }
                  if (rho) {
                    Mhad(NIch(res),NIem("rho"),la1,la2,mu) +=
                      fac*rhofac * (-e/grho_tilde)*F_rho * 
                      ubar(0,la2) *
                      ( iso_s *
                        ( (not sch) ? gamma_null
                          : ( cutoff * vertex3hNrho(g1,g2,g3,spinparity,al,p1+q,mu,-k) *
                              i_ * (gamma_(p1+q)+gamma_unit*mr) * BWs *
                              P3h(p1+q,mr,al,be)*sign_(al)*sign_(be) *
                              vertex3hNpi(g0,spinparity,be,-p1-q,q) ) )
                        + iso_u *
                        ( (not uch) ? gamma_null
                          : ( vertex3hNpi(g0,spinparity,be,p1-k,q) *
                              i_ * (gamma_(p1-k)+gamma_unit*mr) * BWu *
                              P3h(p1-k,mr,be,al)*sign_(al)*sign_(be) *
                              vertex3hNrho(g1,g2,g3,spinparity,al,-p1+k,mu,-k) ) )
                        )
                      * u(0,la1);
                  }
                  if (rho2) {
                    Mhad(NIch(res),NIem("rho"),la1,la2,mu) +=
                      fac*rhofac * (-e/grho_tilde)*(mrho*mrho/M2)*F_rho * 
                      ubar(0,la2) *
                      ( iso_s *
                        ( (not sch) ? gamma_null
                          : ( cutoff * vertex3hNrho(g1,g2,g3,spinparity,al,p1+q,mu,-k) *
                              i_ * (gamma_(p1+q)+gamma_unit*mr) * BWs *
                              P3h(p1+q,mr,al,be)*sign_(al)*sign_(be) *
                              vertex3hNpi(g0,spinparity,be,-p1-q,q) ) )
                        + iso_u *
                        ( (not uch) ? gamma_null
                          : ( vertex3hNpi(g0,spinparity,be,p1-k,q) *
                              i_ * (gamma_(p1-k)+gamma_unit*mr) * BWu *
                              P3h(p1-k,mr,be,al)*sign_(al)*sign_(be) *
                              vertex3hNrho(g1,g2,g3,spinparity,al,-p1+k,mu,-k) ) )
                        )
                      * u(0,la1);
                  }
                }
              }
            }
            //------------  spin-5/2  ------------------------------
            if (spin == 5*half) {
              for (uint al1(0); al1<4; al1++) {
                for (uint al2(0); al2<4; al2++) {
                  for (uint be1(0); be1<4; be1++) {
                    for (uint be2(0); be2<4; be2++) {
                      if (gamma) {
                        Mhad(NIch(res),NIem("gamma"),la1,la2,mu) += 
                          fac*
                          ubar(0,la2) * 
                          ( iso_ds *
                            ( (not sch) ? gamma_null
                              : ( cutoff * vertex5hNgamma(gA,0,0,spinparity,al1,al2,p1+q,mu,-k) *
                                  i_ * (gamma_(p1+q)+gamma_unit*mr) * BWs *
                                  P5h(p1+q,mr,al1,al2,be1,be2)*sign_(al1)*sign_(al2)*sign_(be1)*sign_(be2) * 
                                  vertex5hNpi(g0,spinparity,be1,be2,-p1-q,q) ) )
                            + iso_du *
                            ( (not uch) ? gamma_null
                              : ( vertex5hNpi(g0,spinparity,be1,be2,p1-k,q) *
                                  i_ * (gamma_(p1-k)+gamma_unit*mr) * BWu *
                                  P5h(p1-k,mr,be1,be2,al1,al2)*sign_(al1)*sign_(al2)*sign_(be1)*sign_(be2) * 
                                  vertex5hNgamma(gA,0,0,spinparity,al1,al2,-p1+k,mu,-k) ) )
                            )
                          * u(0,la1);
                      }
                      if (rho) {
                        Mhad(NIch(res),NIem("rho"),la1,la2,mu) += 
                          fac*rhofac * (-e/grho_tilde)*F_rho * 
                          ubar(0,la2) * 
                          ( iso_s *
                            ( (not sch) ? gamma_null
                              : ( cutoff * vertex5hNrho(g1,g2,g3,spinparity,al1,al2,p1+q,mu,-k) *
                                  i_ * (gamma_(p1+q)+gamma_unit*mr) * BWs *
                                  P5h(p1+q,mr,al1,al2,be1,be2)*sign_(al1)*sign_(al2)*sign_(be1)*sign_(be2) * 
                                  vertex5hNpi(g0,spinparity,be1,be2,-p1-q,q) ) )
                            + iso_u *
                            ( (not uch) ? gamma_null
                              : ( vertex5hNpi(g0,spinparity,be1,be2,p1-k,q) *
                                  i_ * (gamma_(p1-k)+gamma_unit*mr) * BWu *
                                  P5h(p1-k,mr,be1,be2,al1,al2)*sign_(al1)*sign_(al2)*sign_(be1)*sign_(be2) * 
                                  vertex5hNrho(g1,g2,g3,spinparity,al1,al2,-p1+k,mu,-k) ) )
                            )
                          * u(0,la1);
                      }
                      if (rho2) {
                        Mhad(NIch(res),NIem("rho"),la1,la2,mu) += 
                          fac*rhofac * (-e/grho_tilde)*(mrho*mrho/M2)*F_rho * 
                          ubar(0,la2) * 
                          ( iso_s *
                            ( (not sch) ? gamma_null
                              : ( cutoff * vertex5hNrho(g1,g2,g3,spinparity,al1,al2,p1+q,mu,-k) *
                                  i_ * (gamma_(p1+q)+gamma_unit*mr) * BWs *
                                  P5h(p1+q,mr,al1,al2,be1,be2)*sign_(al1)*sign_(al2)*sign_(be1)*sign_(be2) * 
                                  vertex5hNpi(g0,spinparity,be1,be2,-p1-q,q) ) )
                            + iso_u *
                            ( (not uch) ? gamma_null
                              : ( vertex5hNpi(g0,spinparity,be1,be2,p1-k,q) *
                                  i_ * (gamma_(p1-k)+gamma_unit*mr) * BWu *
                                  P5h(p1-k,mr,be1,be2,al1,al2)*sign_(al1)*sign_(al2)*sign_(be1)*sign_(be2) * 
                                  vertex5hNrho(g1,g2,g3,spinparity,al1,al2,-p1+k,mu,-k) ) )
                            )
                          * u(0,la1);
                      }
                    }
                  }
                }
              }
            }
          } // if (Config::exists(res))
        } // for (const string& res : resonances)
      } // for (halfint mu(0); mu<4; mu++)
    } // for (halfint la2 : {-half,half})
  } // for (halfint la1 : {-half,half})

  HelicityAmplitudes HA({half,half,_1}, {channels,EMchannels});
  for (halfint la1 : {-half,half}) {
    for (halfint la2 : {-half,half}) {
      for (halfint la : {-_1,_0,_1}) {
        for (string ch : channels) {
          for (string emch : EMchannels) {
            HA({ch,emch},{la1,la2,la}) = 0;
            for (halfint mu(0); mu<4; mu++) {
              HA({ch,emch},{la1,la2,la}) +=
                //Mhad(NIch(ch),NIem(emch),la1,la2,mu) * k(mu) * sign_(mu); // for checking gauge invariance
                Mhad(NIch(ch),NIem(emch),la1,la2,mu) * (epsRe(la)(mu)-i_*epsIm(la)(mu)) * sign_(mu);
               // e/M2 * Mhad(NIch(ch),NIem(emch),la1,la2,mu) * (epsRe(la)(mu)-i_*epsIm(la)(mu)) * sign_(mu);
            }
          }
        }
      }
    }
  }
  //cerr << "check gauge invar (Mhad^mu * eps_mu  -> Mhad^mu * k_mu)  !!!" << endl;
  //cerr << "MUST be changed back !!!" << endl;
  //  cerr << "relative uncertainties have to be set!!!" << endl;
  //  cerr << "return from HelicityAmplitudes" << endl;
  return HA;
}


  /**
     Calculate the production density matrix of the virtual photon
  */
MultiArray<dcomplex> piN_Ngammastar::rho_prod(double costh) const {
  HelicityAmplitudes HA = helicityAmplitudes(costh);
  MultiArray<dcomplex> densityMatrix(idx_s1,idx_s1);
  for (halfint la : {-_1,_0,_1}) {
    for (halfint lap : {-_1,_0,_1}) {
      densityMatrix(la,lap) = 0;
      for (halfint la1 : {-half,half}) {
        for (halfint la2 : {-half,half}) {
          for (string ch : channels) {
            for (string chp : channels) {
              for (string emch : EMchannels) {
                for (string emchp : EMchannels) {
                  densityMatrix(la,lap) += HA({ch,emch},{la1,la2,la}) * conj(HA({chp,emchp},{la1,la2,lap}));
                }
              }
            }
          }
        }
      }
    }
  }

  return densityMatrix;
}

vector<string> piN_Ngammastar::getResonances() const {
  return resonances;
}

vector<string> piN_Ngammastar::getChannels() const {
  return channels;
}

vector<string> piN_Ngammastar::getEMchannels() const {
  return EMchannels;
}
