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

    double util_couple(double Cw_priv, double Cm_priv, double C_pub,int iP, int iL,par_struct* par){
        double love = par->grid_love[iL];
        double power = par->grid_power[iP];

        double Uw = util(Cw_priv,C_pub,woman,par,love);
        double Um = util(Cm_priv,C_pub,man,par,love);

        return power*Uw + (1.0-power)*Um;
    }
}