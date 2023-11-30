#ifndef MAIN
#define PRECOMPUTE
#include "myheader.cpp"
#endif

namespace precompute{

    typedef struct { 
    double power;
    double C_tot;

    par_struct *par;
    } solver_precompute_struct;

    double objfunc_precompute(unsigned n, const double *x, double *grad, void *solver_data_in){
        // unpack
        solver_precompute_struct *solver_data = (solver_precompute_struct *) solver_data_in; 
        
        double C_tot = solver_data->C_tot;
        double power = solver_data->power;
        par_struct *par = solver_data->par;

        double Cw_priv = x[0];  
        double Cm_priv = x[1];
        double C_pub = C_tot - Cw_priv - Cm_priv;

        // weighted utility of choice
        double love = 0.0; // does not matter for optimal allocation
        double val = power*utils::util(Cw_priv,C_pub,woman,par,love) + (1.0-power)*utils::util(Cm_priv,C_pub,man,par,love);

        // return negative of value
        return - val;
    }

    void solve_intraperiod_couple(double* Cw_priv,double* Cm_priv,double* C_pub , double C_tot,double power,par_struct *par){
        
        // setup numerical solver
        solver_precompute_struct* solver_data = new solver_precompute_struct;  
                
        int dim = 2;
        double lb[2],ub[2],x[2];   
        
        auto opt = nlopt_create(NLOPT_LN_BOBYQA, dim); // NLOPT_LD_MMA NLOPT_LD_LBFGS NLOPT_GN_ORIG_DIRECT
        double minf=0.0;

        // search over optimal total consumption, C
        // settings
        solver_data->C_tot = C_tot;         
        solver_data->power = power;
        solver_data->par = par;
        nlopt_set_min_objective(opt, objfunc_precompute, solver_data);   
        nlopt_set_maxeval(opt, 2000);

        // bounds
        lb[0] = 0.0;                
        lb[1] = 0.0;
        ub[0] = solver_data->C_tot;
        ub[1] = solver_data->C_tot;
        nlopt_set_lower_bounds(opt, lb);
        nlopt_set_upper_bounds(opt, ub);

        // optimize
        x[0] = solver_data->C_tot/3.0;
        x[1] = solver_data->C_tot/3.0;
        nlopt_optimize(opt, x, &minf);          
        nlopt_destroy(opt);                 
        
        // unpack
        Cw_priv[0] = x[0];
        Cm_priv[0] = x[1];
        C_pub[0] = C_tot - Cw_priv[0] - Cm_priv[0];

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


    void precompute_margu_single(int i, int gender, par_struct* par){
        double* grid_marg_u_single = par->grid_marg_u_single_w;
        double* grid_marg_u_single_for_inv = par->grid_marg_u_single_w_for_inv;
        if (gender==man){
            grid_marg_u_single = par->grid_marg_u_single_m;
            grid_marg_u_single_for_inv = par->grid_marg_u_single_m_for_inv;
        }

        // calculate marginal utility and inverse marginal utility for EGM
        // utility at current allocation

        if (par->analytic_marg_u_single){
            grid_marg_u_single[i] = utils::marg_util_C(par->grid_C_for_marg_u[i], gender, par);

            // inverse marginal utility: flip the grid of marginal util (such that ascending) and store as new "x-axis" grid
            grid_marg_u_single_for_inv[par->num_marg_u-1 -i] = grid_marg_u_single[i];
        }
        else {
            double delta = 0.0001;
            double util_delta = utils::util_C_single(par->grid_C_for_marg_u[i] + delta,gender,par);
            double util = utils::util_C_single(par->grid_C_for_marg_u[i],gender,par);
            grid_marg_u_single[i] = (util_delta - util)/(delta);
            grid_marg_u_single_for_inv[par->num_marg_u-1 -  i] = grid_marg_u_single[i];
        } //finite difference        
    }


    double util_C_couple(double C_tot, int iP, int iL, par_struct* par){
        double love = par->grid_love[iL];
        double power = par->grid_power[iP];
        double Cw_priv {};
        double Cm_priv {};
        double C_pub {};

        solve_intraperiod_couple(&Cw_priv, &Cm_priv, &C_pub , C_tot,power,par);

        return utils::util_couple(Cw_priv,Cm_priv,C_pub,iP,iL,par);
    }

    void precompute_margu_couple(int i, int iP, par_struct *par, sol_struct *sol){
        double C_tot = par->grid_C_for_marg_u[i];
        int idx = index::index2(iP,i,par->num_power,par->num_marg_u);   

        // calculate marginal utility and inverse marginal utility for EGM
        int iL = 0; // does not matter for the marginal utility     

        // utility at current allocation 
        double util = util_C_couple(C_tot,iP,iL,par);
        par->grid_util[idx] = util; //<-------------- AMO: delete if unnecessary

        // marginal utility
        double delta = 0.0001;
        double util_delta = util_C_couple(C_tot + delta,iP,iL,par);
        par->grid_marg_u[idx] = (util_delta- util)/(delta);

        int idx_flip = index::index2(iP,par->num_marg_u-1 - i,par->num_power,par->num_marg_u);
        par->grid_marg_u_for_inv[idx_flip] = par->grid_marg_u[idx];
    } // precompute func


    void precompute(sol_struct* sol, par_struct* par){
        #pragma omp parallel num_threads(par->threads)      
        {   
            // pre-compute optimal allocation for couple
            #pragma omp for         
            for (int i=0; i<par->num_Ctot; i++){  
                for (int iP=0; iP < par->num_power; iP++){
                    double C_tot = par->grid_Ctot[i];
                    int idx = index::index2(iP,i,par->num_power,par->num_Ctot);         
                    solve_intraperiod_couple(&sol->pre_Ctot_Cw_priv[idx], &sol->pre_Ctot_Cm_priv[idx], &sol->pre_Ctot_C_pub[idx] , C_tot,par->grid_power[iP],par);
                } // power
            } // Ctot
        
            // pre-compute marginal utilities for EGM
            if (par->do_egm){
                #pragma omp for
                for (int i=0; i<par->num_marg_u; i++){  
                    precompute_margu_single(i, woman, par);
                    precompute_margu_single(i, man, par);
                    for (int iP=0; iP < par->num_power; iP++){
                        precompute_margu_couple(i, iP, par, sol);
                    } // power
                } //Ctot
            } // do_egm
        } // parallel
    } // precompute func
}

