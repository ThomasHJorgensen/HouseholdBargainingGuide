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

    double util_couple(double Cw_priv, double Cm_priv, double C_pub,double power, int iL,par_struct* par){
        double love = par->grid_love[iL];

        double Uw = util(Cw_priv,C_pub,woman,par,love);
        double Um = util(Cm_priv,C_pub,man,par,love);

        return power*Uw + (1.0-power)*Um;
    }


    double cons_priv_single(double C_tot,int gender,par_struct *par){
        // closed form solution for intra-period problem of single.
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
        
        return C_tot/(1.0 + pow(alpha2/alpha1,1.0/(1.0-phi) ));
    }

    double util_C_single(double C_tot, int gender, par_struct* par){
        double love = 0.0;
        
        // flow-utility
        double C_priv = cons_priv_single(C_tot,gender,par);
        double C_pub = C_tot - C_priv;
        
        return util(C_priv,C_pub,gender,par,love);
    }

    double marg_util_C(double C_tot, int gender, par_struct* par){
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
        
        double share = 1.0/(1.0 + pow(alpha2/alpha1,1.0/(1.0-phi) ));
        double constant = alpha1*pow(share,phi) + alpha2*pow(1.0-share,phi);
        return phi * pow(C_tot,(1.0-rho)*phi -1.0 ) * pow(constant,1.0 - rho);
    }

    double inv_marg_util_C(double marg_U, int gender, par_struct* par){
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
        
        double share = 1.0/(1.0 + pow(alpha2/alpha1,1.0/(1.0-phi) ));
        double constant = alpha1*pow(share,phi) + alpha2*pow(1.0-share,phi);
        return pow(marg_U / (phi * pow(constant,(1.0-rho))), 1 / ((1-rho)*phi - 1.0));
    }
}