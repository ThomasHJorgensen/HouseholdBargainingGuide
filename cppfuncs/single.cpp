#ifndef MAIN
#define SINGLE
#include "myheader.cpp"
#endif

namespace single {
    typedef struct {
        
        double M;             
        double *V_next;      
        int gender;          
        par_struct *par;      

    } solver_single_struct;

    void intraperiod_allocation(double* C_priv, double* C_pub , double C_tot, int gender,par_struct *par){
        C_priv[0] = utils::cons_priv_single(C_tot,gender,par);
        C_pub[0] = C_tot - C_priv[0];
    }

    // }
    double resources(double A,int gender,par_struct* par) {
        double income = par->inc_w;
        if (gender == man) {
            income = par->inc_m;
        }
        return par->R*A + income;
    }

    double value_of_choice(double C_tot,double M, int gender, double* V_next, par_struct* par){

        // flow-utility
        double Util = utils::util_C_single(C_tot,gender,par);
        
        // continuation value
        double *grid_A = par->grid_Aw; 
        if (gender==man){
            grid_A = par->grid_Am;
        }
        double A = M - C_tot;

        double Vnext = tools::interp_1d(grid_A,par->num_A,V_next,A);
        
        // return discounted sum
        return Util + par->beta*Vnext;
    }

    double objfunc_single(unsigned n, const double *x, double *grad, void *solver_data_in){
        double love = 0.0;

        // unpack
        solver_single_struct *solver_data = (solver_single_struct *) solver_data_in;  
        
        double C_tot = x[0];
        int gender = solver_data->gender;
        double M = solver_data->M;
        par_struct *par = solver_data->par;

        return - value_of_choice(C_tot,M,gender,solver_data->V_next,par);

    }

    void EGM_step_single(int t, sol_struct* sol, par_struct* par){
        // set index for next period (first asset point)
        int idx_next = index::index2(t+1,0,par->T,par->num_A);

        // allocate 
        

        // EGM
        for (int iA_pd=0; iA_pd<par->num_A_pd;iA_pd++){

            // Unpack
            double A_next = par->grid_A_pd[iA_pd];

            // Calculate expected marginal utility
            if (par->do_egm_normal==1){
                // find optimal consumption in next period by interpolation
                double C_next_w = tools::interp_1d(par->grid_Aw, par->num_A,&sol->Cw_tot_single[idx_next], A_next);
                double C_next_m = tools::interp_1d(par->grid_Am, par->num_A,&sol->Cm_tot_single[idx_next], A_next);
                
                // marginal utility of optimal consumption
                sol->EmargUw_single_pd[iA_pd] = par->beta * par->R * utils::marg_util_C(C_next_w, woman, par);
                sol->EmargUm_single_pd[iA_pd] = par->beta * par->R * utils::marg_util_C(C_next_m, man,   par);
            } else{
                // interpolate next period marginal value/consumption
                sol->EmargUw_single_pd[iA_pd] = tools::interp_1d(par->grid_Aw,par->num_A,&sol->marg_Vw_single[idx_next],A_next);
                sol->EmargUm_single_pd[iA_pd] = tools::interp_1d(par->grid_Am,par->num_A,&sol->marg_Vm_single[idx_next],A_next);
            }

            // invert marginal utility by interpolation from pre-computed grid
            sol->C_totw_single_pd[iA_pd] = tools::interp_1d(par->grid_marg_u_single_w_for_inv, par->num_Ctot, par->grid_inv_marg_u, sol->EmargUw_single_pd[iA_pd]);
            sol->C_totm_single_pd[iA_pd] = tools::interp_1d(par->grid_marg_u_single_m_for_inv, par->num_Ctot, par->grid_inv_marg_u, sol->EmargUm_single_pd[iA_pd]);
            
            // endogenous grid over resources
            sol->Mw_single_pd[iA_pd] = sol->C_totw_single_pd[iA_pd] + A_next;
            sol->Mm_single_pd[iA_pd] = sol->C_totm_single_pd[iA_pd] + A_next;
        }

        // interpolate to common grid
        for (int iA=0; iA<par->num_A;iA++){
            int idx = index::index2(t,iA,par->T,par->num_A);

            // create alias for consumption
            double &Cw_tot = sol->Cw_tot_single[idx];
            double &Cm_tot = sol->Cm_tot_single[idx];

            // resources
            double Mw_now = resources(par->grid_Aw[iA], woman, par);
            double Mm_now = resources(par->grid_Am[iA], man, par);

            // total consumption
            Cw_tot = tools::interp_1d(sol->Mw_single_pd, par->num_A_pd, sol->C_totw_single_pd, Mw_now);
            Cm_tot = tools::interp_1d(sol->Mm_single_pd, par->num_A_pd, sol->C_totm_single_pd, Mm_now);

            // woman: if credit constrained
            if (Cw_tot > Mw_now) {
                Cw_tot = Mw_now; // consume all resources

                // marginal value of constrained consumption
                sol->marg_Vw_single[idx] = par->beta*par->R*utils::marg_util_C(Cw_tot, woman, par);
            }
            else{ // if not credit constrained
                // marginal value of unconstrained consumption
                sol->marg_Vw_single[idx] = par->beta*par->R*tools::interp_1d(sol->Mw_single_pd, par->num_A_pd, sol->EmargUw_single_pd, Mw_now);
            }

            // man: if credit constrained
            if (Cm_tot > Mm_now) {
                Cm_tot = Mm_now; // consume all resources

                // marginal value of constrained consumption
                sol->marg_Vm_single[idx] = par->beta*par->R*utils::marg_util_C(Cw_tot, man, par);
            }
            else{ // if not credit constrained
                // marginal value of unconstrained consumption
                sol->marg_Vm_single[idx] = par->beta*par->R*tools::interp_1d(sol->Mm_single_pd, par->num_A_pd, sol->EmargUm_single_pd, Mm_now);
            }

            // allcoate consumption
            intraperiod_allocation(&sol->Cw_priv_single[idx], &sol->Cw_pub_single[idx], Cw_tot, woman, par);
            intraperiod_allocation(&sol->Cm_priv_single[idx], &sol->Cm_pub_single[idx], Cm_tot, man, par);
            
            // calculate values (not used in EGM)
            sol->Vw_single[idx] = value_of_choice(Cw_tot, Mw_now, woman, &sol->Vw_single[idx_next], par);
            sol->Vm_single[idx] = value_of_choice(Cm_tot, Mm_now,   man, &sol->Vm_single[idx_next], par);
        }
    }



    void solve_single(int t,sol_struct *sol,par_struct *par){
        double love = 0.0; // no love for singles 

        // terminal period
        if (t == (par->T-1)){
            for (int iA=0; iA<par->num_A;iA++){
                int idx = index::index2(t,iA,par->T,par->num_A);

                double Aw = par->grid_Aw[iA];
                double Am = par->grid_Am[iA];

                sol->Cw_tot_single[idx] = resources(Aw,woman,par); 
                sol->Cm_tot_single[idx] = resources(Am,man,par); 
                
                intraperiod_allocation(&sol->Cw_priv_single[idx],&sol->Cw_pub_single[idx],sol->Cm_tot_single[idx],woman,par);
                sol->Vw_single[idx] = utils::util(sol->Cw_priv_single[idx],sol->Cw_pub_single[idx],woman,par,love);
                
                intraperiod_allocation(&sol->Cm_priv_single[idx],&sol->Cm_pub_single[idx],sol->Cm_tot_single[idx],man,par);
                sol->Vm_single[idx] = utils::util(sol->Cm_priv_single[idx],sol->Cm_pub_single[idx],man,par,love);

                if (par->do_egm) {
                    sol->marg_Vw_single[idx] = par->beta*par->R*utils::marg_util_C(sol->Cm_tot_single[idx], woman, par);
                    sol->marg_Vm_single[idx] = par->beta*par->R*utils::marg_util_C(sol->Cm_tot_single[idx], man, par);
                }
            }
        } else {
            if (par->do_egm) {
                EGM_step_single(t, sol, par);
            }
            else {
                #pragma omp parallel num_threads(par->threads)
                {

                    // 1. allocate objects for solver
                    solver_single_struct* solver_data = new solver_single_struct;
                    
                    int dim = 1;
                    double lb[1],ub[1],x[1];
                    
                    auto opt = nlopt_create(NLOPT_LN_BOBYQA, dim); // NLOPT_LD_MMA NLOPT_LD_LBFGS NLOPT_GN_ORIG_DIRECT
                    double minf=0.0;

                    // 2. loop over assets
                    #pragma omp for
                    for (int iA=0; iA<par->num_A;iA++){
                        int idx = index::index2(t,iA,par->T,par->num_A);
                        
                        // resources
                        double Aw = par->grid_Aw[iA];
                        double Am = par->grid_Am[iA];
                        
                        double Mw = resources(Aw,woman,par); 
                        double Mm = resources(Am,man,par); 
                        
                        // search over optimal total consumption, C
                        // WOMEN
                        // settings
                        solver_data->M = Mw;
                        solver_data->V_next = &sol->Vw_single[index::index2(t+1,0,par->T,par->num_A)];
                        solver_data->gender = woman;
                        solver_data->par = par;
                        nlopt_set_min_objective(opt, objfunc_single, solver_data); 

                        // bounds
                        lb[0] = 1.0e-8;
                        ub[0] = solver_data->M;
                        nlopt_set_lower_bounds(opt, lb);
                        nlopt_set_upper_bounds(opt, ub);

                        // optimize
                        x[0] = solver_data->M/2.0; 
                        nlopt_optimize(opt, x, &minf); 

                        // store results
                        double Cw = x[0];
                        intraperiod_allocation(&sol->Cw_priv_single[idx],&sol->Cw_pub_single[idx],Cw,woman,par);
                        sol->Vw_single[idx] = -minf;

                        // MEN
                        // settings
                        solver_data->M = Mm;
                        solver_data->V_next = &sol->Vm_single[index::index2(t+1,0,par->T,par->num_A)];
                        solver_data->gender = man;
                        solver_data->par = par;
                        nlopt_set_min_objective(opt, objfunc_single, solver_data);
                            
                        // bounds
                        lb[0] = 1.0e-8;
                        ub[0] = solver_data->M;
                        nlopt_set_lower_bounds(opt, lb);
                        nlopt_set_upper_bounds(opt, ub);

                        // optimize
                        x[0] = solver_data->M/2.0;
                        nlopt_optimize(opt, x, &minf);

                        double Cm = x[0];
                        intraperiod_allocation(&sol->Cm_priv_single[idx],&sol->Cm_pub_single[idx],Cm,man,par);
                        sol->Vm_single[idx] = -minf;            
                        
                    } // iA

                    // 4. destroy optimizer
                    nlopt_destroy(opt);

                } // pragma
            }
        }   
        
    }


}