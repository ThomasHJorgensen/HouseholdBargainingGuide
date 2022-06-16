// functions related to utility and environment.
#ifndef MAIN
#define UTILS
#include "myheader.cpp"
#endif

namespace utils {
    double util(double c_priv,double c_pub,int gender,par_struct *par,double love){
        double rho = par->rho_w;
        double phi = par->phi_w;
        double alpha1 = par->alpha1_w;
        double alpha2 = par->alpha2_w;
        if (gender == man) {
            rho = par->rho_m;
            phi = par->phi_m;
            alpha1 = par->alpha1_m;
            alpha2 = par->alpha2_m;
        }
        
        return pow(alpha1*pow(c_priv,phi) + alpha2*pow(c_pub,phi),1.0-rho)/(1.0-rho) + love;
    }
}