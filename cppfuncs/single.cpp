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

    void EGM_single(int t, int gender, sol_struct* sol, par_struct* par){
        // 1. Setup
        /// a. unpack
        bool const &analytic_inv_marg_u_single {par->analytic_inv_marg_u_single};
        double* const &grid_inv_marg_u {par->grid_inv_marg_u};

        /// b. gender specific variables
        //// o. woman
        double* grid_A {par->grid_Aw};
        double* grid_A_pd {par->grid_Aw_pd};
        double* grid_marg_u_single_for_inv {par->grid_marg_u_single_w_for_inv};
        double* V {sol->Vw_single};
        double* marg_V {sol->marg_Vw_single};
        double* C_tot {sol->Cw_tot_single};
        double* C_priv {sol->Cw_priv_single};
        double* C_pub {sol->Cw_pub_single};
        //// oo. man
        if (gender == man){
            grid_A = par->grid_Am;
            grid_A_pd = par->grid_Am_pd;
            grid_marg_u_single_for_inv = par->grid_marg_u_single_m_for_inv;
            V = sol->Vm_single;
            marg_V = sol->marg_Vm_single;
            C_tot = sol->Cm_tot_single;
            C_priv = sol->Cm_priv_single;
            C_pub = sol->Cm_pub_single;
        }

        /// c. Allocate memory
        double* EmargU_pd {new double[par->num_A_pd]};
        double* C_tot_pd {new double[par->num_A_pd]};
        double* M_pd {new double[par->num_A_pd]};

        // 2. EGM step
        /// setup
        int idx_next = index::single(t+1,0, par);
        int min_point_A = 0;

        for (int iA_pd=0; iA_pd<par->num_A_pd; iA_pd++){

            /// a. get next period assets
            double A_next = grid_A_pd[iA_pd];

            /// b. calculate expected marginal utility
            min_point_A = tools::binary_search(min_point_A, par->num_A, grid_A, A_next);
            EmargU_pd[iA_pd] = tools::interp_1d_index(grid_A, par->num_A, &marg_V[idx_next],A_next, min_point_A);

            /// c. invert marginal utility by interpolation from pre-computed grid
            if (analytic_inv_marg_u_single == 1){
                C_tot_pd[iA_pd] = utils::inv_marg_util_C(EmargU_pd[iA_pd], gender, par);
            } else {
                C_tot_pd[iA_pd] = tools::interp_1d(grid_marg_u_single_for_inv, par->num_marg_u, grid_inv_marg_u, EmargU_pd[iA_pd]);
            }
            /// d. endogenous grid over resources
            M_pd[iA_pd] = C_tot_pd[iA_pd] + A_next;
        }

        // 3. interpolate to common grid
        ///setup
        min_point_A = 0;

        for (int iA=0; iA<par->num_A; iA++){
            int idx = index::single(t,iA, par);

            /// a. calculate resources
            double M_now = resources(grid_A[iA], gender, par);

            /// b. find total consumption
            min_point_A = tools::binary_search(min_point_A,par->num_A_pd, M_pd, M_now);
            C_tot[idx] = tools::interp_1d_index(M_pd, par->num_A_pd, C_tot_pd, M_now, min_point_A);

            /// c. handle credit constraint 
            //// if credit constrained
            if (M_now <= M_pd[0]) {
                ///// o. consume all resources
                C_tot[idx] = M_now; 

                ///// oo. calculate marginal value of constrained consumption
                if (par->analytic_marg_u_single){
                    marg_V[idx] = par->beta * par->R * utils::marg_util_C(C_tot[idx], gender, par);
                }
                else{
                    marg_V[idx] = par->beta * par->R * tools::interp_1d(par->grid_C_for_marg_u, par->num_marg_u, par->grid_marg_u, C_tot[idx]);
                }
            }
            //// if not credit constrained
            else{
                // o. calculate marginal value of unconstrained consumption
                marg_V[idx] = par->beta * par->R * tools::interp_1d_index(M_pd, par->num_A_pd, EmargU_pd, M_now, min_point_A);
            }

            /// d. calculate private and public consumption
            intraperiod_allocation(&C_priv[idx], &C_pub[idx], C_tot[idx], gender, par);
            
            /// e. calculate values (not used in EGM step)
            V[idx] = value_of_choice(C_tot[idx], M_now, gender, &V[idx_next], par);

        }

        // 4. clean up
        delete[] EmargU_pd;
        delete[] C_tot_pd;
        delete[] M_pd;
    }



    void solve_remain_single(int t,sol_struct *sol,par_struct *par){
        double love = 0.0; // no love for singles 

        // terminal period
        if (t == (par->T-1)){
            for (int iA=0; iA<par->num_A;iA++){
                int idx = index::single(t,iA,par);

                double Aw = par->grid_Aw[iA];
                double Am = par->grid_Am[iA];

                sol->Cw_tot_single[idx] = resources(Aw,woman,par); 
                sol->Cm_tot_single[idx] = resources(Am,man,par); 
                
                intraperiod_allocation(&sol->Cw_priv_single[idx],&sol->Cw_pub_single[idx],sol->Cw_tot_single[idx],woman,par);
                sol->Vw_single[idx] = utils::util(sol->Cw_priv_single[idx],sol->Cw_pub_single[idx],woman,par,love);
                
                intraperiod_allocation(&sol->Cm_priv_single[idx],&sol->Cm_pub_single[idx],sol->Cm_tot_single[idx],man,par);
                sol->Vm_single[idx] = utils::util(sol->Cm_priv_single[idx],sol->Cm_pub_single[idx],man,par,love);

                if (par->do_egm) {
                    sol->marg_Vw_single[idx] = par->beta*par->R*utils::marg_util_C(sol->Cw_tot_single[idx], woman, par);
                    sol->marg_Vm_single[idx] = par->beta*par->R*utils::marg_util_C(sol->Cm_tot_single[idx], man, par);
                }
            }
        } else {
            if (par->do_egm) {
                EGM_single(t, woman, sol, par);
                EGM_single(t, man, sol, par);
            }
            else {
                #pragma omp parallel num_threads(par->threads)
                {

                    // 1. allocate objects for solver
                    solver_single_struct* solver_data = new solver_single_struct;
                    
                    int const dim = 1;
                    double lb[dim],ub[dim],x[dim];
                    
                    auto opt = nlopt_create(NLOPT_LN_BOBYQA, dim); // NLOPT_LD_MMA NLOPT_LD_LBFGS NLOPT_GN_ORIG_DIRECT
                    double minf=0.0;

                    // 2. loop over assets
                    #pragma omp for
                    for (int iA=0; iA<par->num_A;iA++){
                        int idx = index::single(t,iA,par);
                        
                        // resources
                        double Aw = par->grid_Aw[iA];
                        double Am = par->grid_Am[iA];
                        
                        double Mw = resources(Aw,woman,par); 
                        double Mm = resources(Am,man,par); 
                        
                        // search over optimal total consumption, C
                        // WOMEN
                        // settings
                        solver_data->M = Mw;
                        solver_data->V_next = &sol->Vw_single[index::single(t+1,0,par)]; // sol->EVw_start_single;
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
                        sol->Cw_tot_single[idx] = x[0];
                        intraperiod_allocation(&sol->Cw_priv_single[idx],&sol->Cw_pub_single[idx],sol->Cw_tot_single[idx],woman,par);
                        sol->Vw_single[idx] = -minf;

                        // MEN
                        // settings
                        solver_data->M = Mm;
                        solver_data->V_next = &sol->Vm_single[index::single(t+1,0,par)]; // sol->EVm_start_single;
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

                        sol->Cm_tot_single[idx] = x[0];
                        intraperiod_allocation(&sol->Cm_priv_single[idx],&sol->Cm_pub_single[idx],sol->Cm_tot_single[idx],man,par);
                        sol->Vm_single[idx] = -minf;       
                        
                    } // iA

                    // 4. destroy optimizer
                    nlopt_destroy(opt);

                } // pragma
            }
        }   
        
    }

    void solve_remain_trans_single(int t,sol_struct *sol,par_struct *par){
        // solve for value of remaining single
        solve_remain_single(t,sol,par);

        // add divorce cost to get value of transition into singlehood
        #pragma omp parallel num_threads(par->threads)
        {
            #pragma omp for
            for (int iA=0; iA<par->num_A;iA++){
                int idx = index::single(t,iA,par);

                sol->Vw_trans_single[idx] = sol->Vw_single[idx] - par->div_cost;
                sol->Vm_trans_single[idx] = sol->Vm_single[idx] - par->div_cost;
                sol->Cw_priv_trans_single[idx] = sol->Cw_priv_single[idx];
                sol->Cm_priv_trans_single[idx] = sol->Cm_priv_single[idx];
                sol->Cw_pub_trans_single[idx] = sol->Cw_pub_single[idx]; 

                if (par->do_egm) {
                    sol->marg_Vw_trans_single[idx] = sol->marg_Vw_single[idx];
                    sol->marg_Vm_trans_single[idx] = sol->marg_Vm_single[idx];
                }
            }
        }
        
    }
    void expected_value_start_single(int t,sol_struct* sol,par_struct* par){
        #pragma omp parallel num_threads(par->threads)
        {
            index::index_couple_struct* idx_couple = new index::index_couple_struct;

            // loop over states
            #pragma omp for
            for (int iA=0; iA<par->num_A;iA++){
                int idx = index::single(t,iA,par);
                // sol->EVw_start_single[idx] = sol->Vw_remain_single[idx]; // not there yet
                // sol->EVm_start_single[idx] = sol->Vm_remain_single[idx];
                // add marginal values

            }
        }
    }

    void expected_value_start_single_new(int t,sol_struct* sol,par_struct* par){
        #pragma omp parallel num_threads(par->threads)
        {
            index::index_couple_struct* idx_couple = new index::index_couple_struct;

            // loop over states
            #pragma omp for
            for (int iA=0; iA<par->num_A;iA++){

                // value of remaining single
                int idx_single = index::single(t,iA,par);
                double V_remain_w = sol->Vw_single[idx_single]; //remain
                double V_remain_m = sol->Vm_single[idx_single];

                // loop over potential partners conditional on meeting a partner
                // love is on the grid, so no need to interpolate in that direction. For wealth we need to.
                double Ev_cond_w = 0.0;
                double Ev_cond_m = 0.0;
                double val_w = 0.0;
                double val_m = 0.0;
                for(int i_love=0;i_love<par->num_love;i_love++){
                    for(int iAp=0;iAp<par->num_A;iAp++){ // partner's wealth

                        // probability of meeting a specific type of partner
                        int idx_A = index::index2(iA,iAp,par->num_A,par->num_A);
                        double prob_A_w = par->prob_partner_A_w[idx_A]; 
                        double prob_A_m = par->prob_partner_A_m[idx_A]; 
                        double prob_love = par->prob_partner_love[i_love]; 
                        double prob_w = prob_A_w*prob_love;
                        double prob_m = prob_A_m*prob_love;

                        // find value of s->m (interpolate in A_tot) and s->s (on grid in A_j)
                        // bargain over consumption
                        // check if remaining single or joining couple
                        // update expected value
                        

                        // OLD:
                        // find relevant value function 
                        int idx_trans = index::trans_to_couple(t,i_love,iA,par);
                        int power_idx = sol->power_idx_trans[idx_trans];
                        
                        if (power_idx>=0){
                            // TODO: interpolate: note there needs to be done something about wealth! The calculation of the value of transitioning might be move to here!
                            double Aw_tot = par->grid_Aw[iA] + par->grid_Aw[iAp];
                            double Am_tot = par->grid_Am[iA] + par->grid_Am[iAp]; 
                            int idx_interp = index::trans_to_couple(t,i_love,0,par);
                            val_w = tools::interp_1d(par->grid_A,par->num_A,&sol->Vw_trans_couple[idx_interp],Aw_tot);
                            val_m = tools::interp_1d(par->grid_A,par->num_A,&sol->Vm_trans_couple[idx_interp],Am_tot);
                        
                        } else {
                            val_w = V_remain_w;
                            val_m = V_remain_m;
                        }

                        // expected value conditional on meeting a partner
                        Ev_cond_w += prob_w*val_w;
                        Ev_cond_m += prob_m*val_m;

                    }
                }

                // expected value of starting single
                double p_meet = par->prob_repartner[t]; 
                double Ev_w = p_meet*Ev_cond_w + (1.0-p_meet)*V_remain_w;
                double Ev_m = p_meet*Ev_cond_m + (1.0-p_meet)*V_remain_m;

                sol->EVw_start_single[idx_single] = Ev_w;
                sol->EVm_start_single[idx_single] = Ev_m;

                
            }

        } // pragma
    }


}