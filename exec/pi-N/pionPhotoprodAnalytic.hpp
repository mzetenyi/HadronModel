#ifndef PIONPHOTOPRODANALYTIC
#define PIONPHOTOPRODANALYTIC

#include "Vectors.hpp"
using namespace Vectors;
#include "Kinema.hpp"

class pionPhotoprodAnalytic {
 private:
  double srt;
  double mN;
  double mpi;
  Kinema2 KINin;
  Kinema2 KINout;

 public:
  pionPhotoprodAnalytic(double srt);
  ~pionPhotoprodAnalytic();
  double MSQR(double costh);
  double dsig_dcosth(double costh);
  double sigtot();
};

#endif // PIONPHOTOPRODANALYTIC