#include "Isospin.hpp"
#include <cmath>
#include <iostream>


double isospin_1h1h1(int Q1h1, int Q1h2, int Q1) {
    if (Q1==0) {
        if (Q1h1==0) return -1;
        return 1;
    }
    return sqrt(2.);
    //return (Q1==0) ? ((Q1h1==0) ? -1 : 1) : sqrt(2.);
}

double isospin_3h1h1(int Q3h, int Q1h, int Q1) {
    if (Q3h==2) return -1;
    if (Q3h==1) {
        if (Q1==0) return sqrt(2./3.);
        if (Q1==1) return -sqrt(1./3.);
    }
    if (Q3h==0) {
        if (Q1==-1) return sqrt(1./3.);
        if (Q1==0) return sqrt(2./3.);
    }
    if (Q3h==-1) return 1;
    std::cerr << "isospin_3h1h1: undefined for (Q3h,Q1h,Q1) =  (" << Q3h << "," << Q1h << "," << Q1 << ")"<< std::endl;
    exit(0);
}

double isospin_1h3h1(int Q3h, int Q1h, int Q1) {
    if (Q3h==-2) return -1;
    if (Q3h==-1) {
        if (Q1==0) return sqrt(2./3.);
        if (Q1==1) return -sqrt(1./3.);
    }
    if (Q3h==0) {
        if (Q1==-1) return sqrt(1./3.);
        if (Q1==0) return sqrt(2./3.);
    }
    if (Q3h==1) return 1;
    std::cerr << "isospin_1h3h1: undefined for (Q1h,Q3h,Q1) =  (" << Q1h << "," << Q3h << "," << Q1 << ")"<< std::endl;
    exit(0);
}
