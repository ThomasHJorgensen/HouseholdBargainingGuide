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

        // clip and penalty
        double penalty = 0.0;
        if(Cw_priv < 1.0e-8){
            penalty += 1000.0*(Cw_priv*Cw_priv);
            Cw_priv = 1.0e-6;
        }
        if(Cm_priv < 1.0e-8){
            penalty += 1000.0*(Cm_priv*Cm_priv);
            Cm_priv = 1.0e-6;
        }

        double C_pub = C_tot - Cw_priv - Cm_priv;

        // weighted utility of choice
        double love = 0.0; // does not matter for optimal allocation
        double val = power*utils::util(Cw_priv,C_pub,woman,par,love) + (1.0-power)*utils::util(Cm_priv,C_pub,man,par,love);

        // return negative of value
        return - val + penalty;
    }

    void solve_intraperiod_couple(double* Cw_priv,double* Cm_priv,double* C_pub , double C_tot,double power,par_struct *par, double start_Cw_priv, double start_Cm_priv, double ftol = 1.0e-6, double xtol = 1.0e-5){
        // setup numerical solver
        solver_precompute_struct* solver_data = new solver_precompute_struct;  
                
        int const dim = 2;
        double lb[dim],ub[dim],x[dim];   
        
        auto opt = nlopt_create(NLOPT_LN_BOBYQA, dim); // NLOPT_LD_MMA NLOPT_LD_LBFGS NLOPT_GN_ORIG_DIRECT
        double minf=0.0;

        // search over optimal total consumption, C
        // settings
        solver_data->C_tot = C_tot;         
        solver_data->power = power;
        solver_data->par = par;
        nlopt_set_min_objective(opt, objfunc_precompute, solver_data);   
        nlopt_set_maxeval(opt, 2000);
        nlopt_set_ftol_rel(opt, ftol);
        nlopt_set_xtol_rel(opt, xtol);

        // bounds
        lb[0] = 1.0e-6;                
        lb[1] = 1.0e-6;
        ub[0] = solver_data->C_tot;
        ub[1] = solver_data->C_tot;

        if (lb[0] >= ub[0]){
            
            Cw_priv[0] = ub[0]*0.33;
            Cm_priv[0] = ub[1]*0.33;
            C_pub[0] = C_tot - Cw_priv[0] - Cm_priv[0];

        } else {
            nlopt_set_lower_bounds(opt, lb);
            nlopt_set_upper_bounds(opt, ub);

            x[0] = start_Cw_priv;
            x[1] = start_Cm_priv;

            nlopt_optimize(opt, x, &minf);                          
            
            // unpack
            Cw_priv[0] = x[0];
            Cm_priv[0] = x[1];
            C_pub[0] = C_tot - Cw_priv[0] - Cm_priv[0];
        }

        nlopt_destroy(opt); 
        delete solver_data;

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


    void precompute_cons_interp_single(int i, int gender, par_struct* par){
        double* grid_marg_u_single = par->grid_marg_u_single_w;
        double* grid_marg_u_single_for_inv = par->grid_marg_u_single_w_for_inv;
        if (gender==man){
            grid_marg_u_single = par->grid_marg_u_single_m;
            grid_marg_u_single_for_inv = par->grid_marg_u_single_m_for_inv;
        }

        // calculate marginal utility and inverse marginal utility for EGM
        double delta = 0.0001;
        double util_delta = utils::util_C_single(par->grid_C_for_marg_u[i] + delta,gender,par);
        double util = utils::util_C_single(par->grid_C_for_marg_u[i],gender,par);
        grid_marg_u_single[i] = (util_delta - util)/(delta);
        grid_marg_u_single_for_inv[par->num_marg_u-1 -  i] = grid_marg_u_single[i];
    
    }

    void intraperiod_allocation(double* Cw_priv, double* Cm_priv, double* C_pub , double C_tot,int iP,sol_struct *sol,par_struct *par){
        // This function is almost identical to one used in pre-computation... to get "util_C_couple"
        if(par->precompute_intratemporal){
            // interpolate pre-computed solution 
            int idx = index::index2(iP,0,par->num_power,par->num_Ctot); 
            int j1 = tools::binary_search(0,par->num_Ctot,par->grid_Ctot,C_tot);

            Cw_priv[0] = tools::interp_1d_index(par->grid_Ctot,par->num_Ctot,&sol->pre_Ctot_Cw_priv[idx],C_tot,j1);
            Cm_priv[0] = tools::interp_1d_index(par->grid_Ctot,par->num_Ctot,&sol->pre_Ctot_Cm_priv[idx],C_tot,j1);
            C_pub[0] = C_tot - Cw_priv[0] - Cm_priv[0];

        } else {
            // solve intertemporal problem
            double start_Cw_priv = C_tot/3.0; // could be inputs in the pointers
            double start_Cm_priv = C_tot/3.0;
            solve_intraperiod_couple(Cw_priv,Cm_priv,C_pub,C_tot,par->grid_power[iP],par,start_Cw_priv,start_Cm_priv);
        
        }
    }

    void intraperiod_allocation_sim(double* Cw_priv, double* Cm_priv, double* C_pub , double C_tot,double power,sol_struct *sol,par_struct *par){
        if(par->precompute_intratemporal){
            // interpolate pre-computed solution in both power and C_tot, different from solution
            int idx = index::index2(0,0,par->num_power,par->num_Ctot); 
            tools::interp_2d_2out(par->grid_power,par->grid_Ctot,par->num_power,par->num_Ctot,&sol->pre_Ctot_Cw_priv[idx],&sol->pre_Ctot_Cm_priv[idx],power,C_tot,Cw_priv,Cm_priv);

            C_pub[0] = C_tot - Cw_priv[0] - Cm_priv[0];
        } else {
            
            // solve intertemporal problem
            double start_Cw_priv = C_tot/3.0; // could be inputs in the pointers
            double start_Cm_priv = C_tot/3.0;
            precompute::solve_intraperiod_couple(Cw_priv,Cm_priv,C_pub,C_tot,power,par,start_Cw_priv,start_Cm_priv);
        
        }

    }


    EXPORT double util_C_couple(double C_tot, int iP, int iL, par_struct* par, sol_struct* sol, double* Cw_priv, double* Cm_priv, double start_Cw_priv, double start_Cm_priv){
        double love = par->grid_love[iL];
        double power = par->grid_power[iP];
        double C_pub = 0.0;
        intraperiod_allocation(Cw_priv,Cm_priv,&C_pub ,C_tot,iP,sol,par);

        return utils::util_couple(*Cw_priv,*Cm_priv,C_pub,power,iL,par);
    }


    EXPORT double marg_util_C_couple(double C_tot, int iP, par_struct* par, sol_struct* sol, double start_Cw_priv, double start_Cm_priv){
        // baseline utility (could be passed as argument to avoid recomputation of utility at C_tot)
        int iL = 0; // does not matter for the marginal utility   

        double Cw_priv = 0.0; 
        double Cm_priv = 0.0;

        double util = util_C_couple(C_tot,iP,iL, par, sol, &Cw_priv,&Cm_priv, start_Cw_priv, start_Cm_priv); // this will update Cw_priv, Cm_priv

        // forward difference
        double delta = 0.0001;
        start_Cw_priv = Cw_priv; //updated with previous solution 
        start_Cm_priv = Cm_priv; 
        double util_delta = util_C_couple(C_tot + delta,iP,iL,par,sol, &Cw_priv,&Cm_priv, start_Cw_priv, start_Cm_priv); 
        return (util_delta - util)/delta;
    }

    void precompute_cons_interp_couple(int i, int iP, par_struct *par, sol_struct *sol){
        
        double C_tot = par->grid_C_for_marg_u[i];
        int idx = index::index2(iP,i,par->num_power,par->num_marg_u);   

        // calculate marginal utility and inverse marginal utility for EGM
        double start_Cw_priv = C_tot/3.0;
        double start_Cm_priv = C_tot/3.0;

        par->grid_marg_u[idx] = marg_util_C_couple(C_tot,iP,par, sol, start_Cw_priv, start_Cm_priv);

        int idx_flip = index::index2(iP,par->num_marg_u-1 - i,par->num_power,par->num_marg_u);
        par->grid_marg_u_for_inv[idx_flip] = par->grid_marg_u[idx];

    } // precompute func


    EXPORT void precompute(sol_struct* sol, par_struct* par){
        #pragma omp parallel num_threads(par->threads)      
        {   
            // pre-compute optimal allocation for couple
            if(par->precompute_intratemporal){    
                #pragma omp for     
                for (int i=0; i<par->num_Ctot; i++){  
                    double C_tot = par->grid_Ctot[i];
                    for (int iP=0; iP < par->num_power; iP++){
                        int idx = index::index2(iP,i,par->num_power,par->num_Ctot);         

                        double start_Cm_priv = C_tot/3.0;
                        double start_Cw_priv = C_tot/3.0;
                        double ftol = 1.0e-10;
                        double xtol = 1.0e-9;
                        if (iP>1){ //Note: using previous starting values does not work super well for iP=1
                            int idx_minus = index::index2(iP-1,i,par->num_power,par->num_Ctot);
                            start_Cm_priv = sol->pre_Ctot_Cm_priv[idx_minus];
                            start_Cw_priv = sol->pre_Ctot_Cw_priv[idx_minus];
                            ftol = 1.0e-6;
                            xtol = 1.0e-5;
                        }
                        solve_intraperiod_couple(&sol->pre_Ctot_Cw_priv[idx], &sol->pre_Ctot_Cm_priv[idx], &sol->pre_Ctot_C_pub[idx] , C_tot,par->grid_power[iP],par, start_Cw_priv, start_Cm_priv, ftol, xtol);
                    } // power
                } // Ctot
            }
        
            // pre-compute marginal utilities and consumption interpolator for EGM
            if( (par->do_egm) & (strcmp(par->interp_method,"numerical")!=0)){
                #pragma omp for
                for (int i=0; i<par->num_marg_u; i++){ 
                    precompute_cons_interp_single(i, woman, par);
                    precompute_cons_interp_single(i, man, par);
                    for (int iP=0; iP < par->num_power; iP++){
                        precompute_cons_interp_couple(i, iP, par, sol);
                    } // power
                } //Ctot
            } // do_egm
        } // parallel
    } // precompute func

    // numerical inverse marginal utility
    typedef struct { 
        double margU;
        int iP;
        int gender;
        par_struct *par;
        sol_struct *sol;
        bool do_print;

        double guess_Cw_priv;
        double guess_Cm_priv;
    } solver_inv_struct;

    double obj_inv_marg_util_couple(unsigned n, const double *x, double *grad, void *solver_data_in){
         // unpack
        solver_inv_struct *solver_data = (solver_inv_struct *) solver_data_in; 
        
        double C_tot = x[0];
        double margU = solver_data->margU;
        int iP = solver_data->iP;
        bool do_print = solver_data->do_print;
        par_struct *par = solver_data->par;
        sol_struct *sol = solver_data->sol;

        // clip
        double penalty = 0.0;
        if (C_tot <= 0.0) {
            penalty += 1000.0*C_tot*C_tot;
            C_tot = 1.0e-6;
        }

        // return squared difference
        double start_Cw_priv = solver_data->guess_Cw_priv;//C_tot/3.0;
        double start_Cm_priv = solver_data->guess_Cm_priv;//C_tot/3.0;
        double diff = marg_util_C_couple(C_tot,iP,par,sol, start_Cw_priv, start_Cm_priv) - margU;

        if (do_print){
            logs::write("inverse_log.txt",1,"C_tot: %f, diff: %f, penalty: %f\n",C_tot,diff,penalty);
        }
        return diff*diff + penalty;

    }

    EXPORT double inv_marg_util_couple(double margU, int iP,par_struct* par, sol_struct* sol, double guess_Ctot, double guess_Cw_priv, double guess_Cm_priv,bool do_print=false ){
        // setup numerical solver
        solver_inv_struct* solver_data = new solver_inv_struct;  
                
        int const dim = 1;
        double lb[dim],ub[dim],x[dim];   
        
        auto opt = nlopt_create(NLOPT_LN_BOBYQA, dim); // NLOPT_LD_MMA NLOPT_LD_LBFGS NLOPT_GN_ORIG_DIRECT NLOPT_LN_BOBYQA
        double minf=0.0;

        // search over optimal total consumption, C
        // settings
        solver_data->margU = margU;         
        solver_data->iP = iP;
        solver_data->par = par;
        solver_data->sol = sol;
        solver_data->do_print = do_print;    

        solver_data->guess_Cw_priv = guess_Cw_priv;
        solver_data->guess_Cm_priv = guess_Cm_priv;


        if (do_print){
            logs::write("inverse_log.txt",0,"margU: %f\n",margU);
        }

        nlopt_set_min_objective(opt, obj_inv_marg_util_couple, solver_data);   
        nlopt_set_maxeval(opt, 2000);
        nlopt_set_ftol_rel(opt, 1.0e-6);
        nlopt_set_xtol_rel(opt, 1.0e-5);

        // bounds
        lb[0] = 0.0;  
        ub[0] = 2.0*par->max_Ctot;
        nlopt_set_lower_bounds(opt, lb);
        nlopt_set_upper_bounds(opt, ub);

        // optimize
        x[0] = guess_Ctot; 
        nlopt_optimize(opt, x, &minf);          
        nlopt_destroy(opt);   

        delete solver_data;              
        
        // return consumption value
        return x[0];
        
    }

    // Could be deleted, or used to compare with analytical.
    double obj_inv_marg_util_single(unsigned n, const double *x, double *grad, void *solver_data_in){
         // unpack
        solver_inv_struct *solver_data = (solver_inv_struct *) solver_data_in; 
        
        double C_tot = x[0];
        double margU = solver_data->margU;
        int gender = solver_data->gender;
        bool do_print = solver_data->do_print;
        par_struct *par = solver_data->par;

        // clip
        double penalty = 0.0;
        if (C_tot <= 0.0) {
            penalty += 1000.0*C_tot*C_tot;
            C_tot = 1.0e-6;
        }

        // return squared difference (using analytical marginal utility)
        double diff = utils::marg_util_C(C_tot, gender, par) - margU;

        if (do_print){
            logs::write("inverse_log.txt",1,"C_tot: %f, diff: %f, penalty: %f\n",C_tot,diff,penalty);
        }
        return diff*diff + penalty;

    }

    EXPORT double inv_marg_util_single(double margU, int gender, par_struct* par, double guess = 3.0, bool do_print=false){
        // setup numerical solver
        solver_inv_struct* solver_data = new solver_inv_struct;  
                
        int const dim = 1;
        double lb[dim],ub[dim],x[dim];   
        
        auto opt = nlopt_create(NLOPT_LN_BOBYQA, dim); // NLOPT_LD_MMA NLOPT_LD_LBFGS NLOPT_GN_ORIG_DIRECT NLOPT_LN_BOBYQA
        double minf=0.0;

        // search over optimal total consumption, C
        // settings
        solver_data->margU = margU;  
        solver_data->gender = gender;         
        solver_data->par = par;
        solver_data->do_print = do_print;        

        if (do_print){
            logs::write("inverse_log.txt",0,"margU: %f\n",margU);
        }

        nlopt_set_min_objective(opt, obj_inv_marg_util_single, solver_data);   
        nlopt_set_maxeval(opt, 2000);
        nlopt_set_ftol_rel(opt, 1.0e-6);
        nlopt_set_xtol_rel(opt, 1.0e-5);

        // bounds
        lb[0] = 0.0;  
        ub[0] = 2.0*par->max_Ctot;
        nlopt_set_lower_bounds(opt, lb);
        nlopt_set_upper_bounds(opt, ub);

        // optimize
        x[0] = guess; 
        nlopt_optimize(opt, x, &minf);          
        nlopt_destroy(opt);                 
        
        delete solver_data;

        // return consumption value
        return x[0];
        
    }
}

