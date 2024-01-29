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

    double value_of_choice_single_to_single(double C_tot,double M, int gender, double* V_next, par_struct* par){

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

    double objfunc_single_to_single(unsigned n, const double *x, double *grad, void *solver_data_in){
        double love = 0.0;

        // unpack
        solver_single_struct *solver_data = (solver_single_struct *) solver_data_in;  
        
        double C_tot = x[0];
        int gender = solver_data->gender;
        double M = solver_data->M;
        par_struct *par = solver_data->par;

        return - value_of_choice_single_to_single(C_tot,M,gender,solver_data->V_next,par);

    }

    void EGM_single_to_single(int t, int gender, sol_struct* sol, par_struct* par){
        // 1. Setup
        /// a. unpack
        bool const &analytic_inv_marg_u_single {par->analytic_inv_marg_u_single};
        double* const &grid_inv_marg_u {par->grid_inv_marg_u};

        /// b. gender specific variables
        //// o. woman
        double* grid_A {par->grid_Aw};
        double* grid_A_pd {par->grid_Aw_pd};
        double* grid_marg_u_single_for_inv {par->grid_marg_u_single_w_for_inv};
        double* V {sol->Vw_single_to_single};
        double* EV {sol->EVw_start_as_single};
        double* margV {sol->EmargVw_start_as_single};
        double* C_tot {sol->Cw_tot_single_to_single};
        double* C_priv {sol->Cw_priv_single_to_single};
        double* C_pub {sol->Cw_pub_single_to_single};
        //// oo. man
        if (gender == man){
            grid_A = par->grid_Am;
            grid_A_pd = par->grid_Am_pd;
            grid_marg_u_single_for_inv = par->grid_marg_u_single_m_for_inv;
            V = sol->Vm_single_to_single;
            EV = sol->EVm_start_as_single;
            margV = sol->EmargVm_start_as_single;
            C_tot = sol->Cm_tot_single_to_single;
            C_priv = sol->Cm_priv_single_to_single;
            C_pub = sol->Cm_pub_single_to_single;
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
            EmargU_pd[iA_pd] = par->beta*tools::interp_1d_index(grid_A, par->num_A, &margV[idx_next],A_next, min_point_A);

            /// c. invert marginal utility by interpolation from pre-computed grid
            if (analytic_inv_marg_u_single == 1){
                C_tot_pd[iA_pd] = utils::inv_marg_util_C(EmargU_pd[iA_pd], gender, par);
            } else {
                if (par->interp_inverse){
                    C_tot_pd[iA_pd] = 1.0/tools::interp_1d(grid_marg_u_single_for_inv, par->num_marg_u, grid_inv_marg_u, EmargU_pd[iA_pd]); // AMO: invert
                }
                else{
                    C_tot_pd[iA_pd] = tools::interp_1d(grid_marg_u_single_for_inv, par->num_marg_u, grid_inv_marg_u, EmargU_pd[iA_pd]);
                }
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

                // ///// oo. calculate marginal value of constrained consumption
                // if (par->analytic_marg_u_single){
                //     margV[idx] = par->R * utils::marg_util_C(C_tot[idx], gender, par);
                // }
                // else{
                //     margV[idx] = par->R * tools::interp_1d(par->grid_C_for_marg_u, par->num_marg_u, par->grid_marg_u, C_tot[idx]);
                // }
            }
            //// if not credit constrained
            // else{
            //     // o. calculate marginal value of unconstrained consumption
            //     margV[idx] = par->R * tools::interp_1d_index(M_pd, par->num_A_pd, EmargU_pd, M_now, min_point_A);
            // }

            /// d. calculate private and public consumption
            intraperiod_allocation(&C_priv[idx], &C_pub[idx], C_tot[idx], gender, par);
            
            /// e. calculate values (not used in EGM step)
            V[idx] = value_of_choice_single_to_single(C_tot[idx], M_now, gender, &EV[idx_next], par);

        }

        // 4. clean up
        delete[] EmargU_pd;
        delete[] C_tot_pd;
        delete[] M_pd;
    }


    void calc_marginal_value_single(int t, int gender, sol_struct* sol, par_struct* par){

        // unpack
        int const &num_A = par->num_A;

        // set index
        int idx = index::index2(t,0,par->T,par->num_A);

        // gender specific variables
        double* grid_A = par->grid_Aw;
        double* margV = &sol->EmargVw_start_as_single[idx];
        double* V      = &sol->EVw_start_as_single[idx];

        if (gender == man){
            grid_A = par->grid_Am;
            margV = &sol->EmargVm_start_as_single[idx];
            V      = &sol->EVm_start_as_single[idx];
        }

        // approximate marginal value by finite diff
        for (int iA=0; iA<num_A-1; iA++){
            // Setup indices
            int iA_plus = iA + 1;

            // Calculate finite difference
            margV[iA] = V[iA_plus] / (grid_A[iA_plus] - grid_A[iA]) - V[iA] / (grid_A[iA_plus] - grid_A[iA]);

            // Extrapolate gradient in last point
            if (iA == num_A-2){
                margV[iA_plus] = margV[iA];
            }
        }
    }



    void solve_single_to_single(int t,sol_struct *sol,par_struct *par){
        double love = 0.0; // no love for singles 

        // terminal period
        if (t == (par->T-1)){
            for (int iA=0; iA<par->num_A;iA++){
                int idx = index::single(t,iA,par);

                double Aw = par->grid_Aw[iA];
                double Am = par->grid_Am[iA];

                sol->Cw_tot_single_to_single[idx] = resources(Aw,woman,par); 
                sol->Cm_tot_single_to_single[idx] = resources(Am,man,par); 
                
                intraperiod_allocation(&sol->Cw_priv_single_to_single[idx],&sol->Cw_pub_single_to_single[idx],sol->Cw_tot_single_to_single[idx],woman,par);
                sol->Vw_single_to_single[idx] = utils::util(sol->Cw_priv_single_to_single[idx],sol->Cw_pub_single_to_single[idx],woman,par,love);
                
                intraperiod_allocation(&sol->Cm_priv_single_to_single[idx],&sol->Cm_pub_single_to_single[idx],sol->Cm_tot_single_to_single[idx],man,par);
                sol->Vm_single_to_single[idx] = utils::util(sol->Cm_priv_single_to_single[idx],sol->Cm_pub_single_to_single[idx],man,par,love);

            }
        } else {
            if (par->do_egm) {
                EGM_single_to_single(t, woman, sol, par);
                EGM_single_to_single(t, man, sol, par);
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
                        solver_data->V_next = &sol->EVw_start_as_single[index::single(t+1,0,par)]; // sol->EVw_start_single;
                        solver_data->gender = woman;
                        solver_data->par = par;
                        nlopt_set_min_objective(opt, objfunc_single_to_single, solver_data); 

                        // bounds
                        lb[0] = 1.0e-8;
                        ub[0] = solver_data->M;
                        nlopt_set_lower_bounds(opt, lb);
                        nlopt_set_upper_bounds(opt, ub);

                        // optimize
                        x[0] = solver_data->M/2.0; 
                        nlopt_optimize(opt, x, &minf); 

                        // store results
                        sol->Cw_tot_single_to_single[idx] = x[0];
                        intraperiod_allocation(&sol->Cw_priv_single_to_single[idx],&sol->Cw_pub_single_to_single[idx],sol->Cw_tot_single_to_single[idx],woman,par);
                        sol->Vw_single_to_single[idx] = -minf;


                        // MEN
                        // settings
                        solver_data->M = Mm;
                        solver_data->V_next = &sol->EVm_start_as_single[index::single(t+1,0,par)]; // sol->EVm_start_singl;
                        solver_data->gender = man;
                        solver_data->par = par;
                        nlopt_set_min_objective(opt, objfunc_single_to_single, solver_data);
                            
                        // bounds
                        lb[0] = 1.0e-8;
                        ub[0] = solver_data->M;
                        nlopt_set_lower_bounds(opt, lb);
                        nlopt_set_upper_bounds(opt, ub);

                        // optimize
                        x[0] = solver_data->M/2.0;
                        nlopt_optimize(opt, x, &minf);

                        sol->Cm_tot_single_to_single[idx] = x[0];
                        intraperiod_allocation(&sol->Cm_priv_single_to_single[idx],&sol->Cm_pub_single_to_single[idx],sol->Cm_tot_single_to_single[idx],man,par);
                        sol->Vm_single_to_single[idx] = -minf;       
                        
                    } // iA

                    // 4. destroy optimizer
                    nlopt_destroy(opt);

                } // pragma
            }
        }   
        
    }


    void solve_couple_to_single(int t, sol_struct *sol, par_struct *par) {
        #pragma omp parallel num_threads(par->threads)
        {
            #pragma omp for
            for (int iA=0; iA<par->num_A;iA++){
                int idx = index::single(t,iA,par);

                sol->Vw_couple_to_single[idx] = sol->Vw_single_to_single[idx] - par->div_cost;
                sol->Vm_couple_to_single[idx] = sol->Vm_single_to_single[idx] - par->div_cost;
                sol->Cw_priv_couple_to_single[idx] = sol->Cw_priv_single_to_single[idx];
                sol->Cm_priv_couple_to_single[idx] = sol->Cm_priv_single_to_single[idx];
                sol->Cw_pub_couple_to_single[idx] = sol->Cw_pub_single_to_single[idx]; 
                sol->Cm_pub_couple_to_single[idx] = sol->Cm_pub_single_to_single[idx]; 
                sol->Cw_tot_couple_to_single[idx] = sol->Cw_tot_single_to_single[idx]; 
                sol->Cm_tot_couple_to_single[idx] = sol->Cm_tot_single_to_single[idx]; 
            }
        }
    }

 int calc_initial_bargaining_weight(int t, int iL, int iAw, int iAm, sol_struct *sol, par_struct *par){
        // a. value of being single
        int idx_single_w = index::single(t,iAw,par);
        int idx_single_m = index::single(t,iAm,par);

        double Vw_single = sol->Vw_single_to_single[idx_single_w];
        double Vm_single = sol->Vm_single_to_single[idx_single_m];

        // b. Setup values for being in couple
        double Vw_single_to_couple = 0.0;
        double Vm_single_to_couple = 0.0;
        double nash_surplus = 0.0;

        int max_idx = -1;
        double max_nash_surplus = 0.0; 
        double A_tot = par->grid_Aw[iAw] + par->grid_Am[iAm];

        double Sw = 0;
        double Sm = 0;

        int iA = tools::binary_search(0, par->num_A, par->grid_A, A_tot);

        // c. loop over bargaining weights
        for (int iP=0; iP < par->num_power; iP++){
            int idx_interp = index::couple(t, iP, iL, 0, par);;
            Vw_single_to_couple = tools::interp_1d_index(par->grid_A, par->num_A, &sol->Vw_single_to_couple[idx_interp], A_tot, iA);
            Vm_single_to_couple = tools::interp_1d_index(par->grid_A, par->num_A, &sol->Vm_single_to_couple[idx_interp], A_tot, iA);
            Sw = Vw_single_to_couple - Vw_single;
            Sm = Vm_single_to_couple - Vm_single;

            // c.1. find power idx that maxes Nash surplus
            if ((Sw>0) & (Sm>0)){
                nash_surplus = Sw*Sm;
                if (nash_surplus > max_nash_surplus){
                    max_nash_surplus = nash_surplus;
                    max_idx = iP;
                }
            }
        }
        return max_idx;
    }

    void expected_value_start_single(int t,sol_struct* sol,par_struct* par){
        #pragma omp parallel num_threads(par->threads)
        {
            index::index_couple_struct* idx_couple = new index::index_couple_struct;

            // a. calculate initial bargaining weights
            // loop over states
            #pragma omp for
            for (int iAw=0; iAw<par->num_A;iAw++){
                for (int iAm=0; iAm<par->num_A;iAm++){
                    // only calculate if match has positive probability of happening
                    int idx_w = index::index2(iAw,iAm,par->num_A,par->num_A);
                    double pw = par->prob_partner_A_w[idx_w]; // woman's prob of meeting man
                    int idx_m = index::index2(iAm,iAw,par->num_A,par->num_A);
                    double pm = par->prob_partner_A_m[idx_w]; // man's prob of meeting

                    if ((pw > 0.0) | (pm > 0.0)) {
                        for (int iL=0; iL<par->num_love;iL++){
                            int idx = index::index4(t,iL,iAw,iAm,par->T,par->num_love,par->num_A,par->num_A);
                            sol->initial_power_idx[idx] = calc_initial_bargaining_weight(t, iL, iAw, iAm, sol, par);
                        } // love
                    } // if
                } // iAm
            } // iAw

            // b. Loop over states
            #pragma omp for
            for (int iA=0; iA<par->num_A;iA++){

                // value of remaining single
                int idx_single = index::single(t,iA,par);
                double Vw_single_to_single = sol->Vw_single_to_single[idx_single]; //remain
                double Vm_single_to_single = sol->Vm_single_to_single[idx_single];

                // b.1. loop over potential partners conditional on meeting a partner
                double Ev_cond_w = 0.0;
                double Ev_cond_m = 0.0;
                double val_w = 0.0;
                double val_m = 0.0;
                for(int iL=0;iL<par->num_love;iL++){
                    for(int iAp=0;iAp<par->num_A;iAp++){ // partner's wealth 

                        // b.1.1. probability of meeting a specific type of partner
                        int idx_A = index::index2(iA,iAp,par->num_A,par->num_A);
                        double prob_A_w = par->prob_partner_A_w[idx_A]; 
                        double prob_A_m = par->prob_partner_A_m[idx_A]; 
                        double prob_love = par->prob_partner_love[iL]; 
                        double prob_w = prob_A_w*prob_love;
                        double prob_m = prob_A_m*prob_love;

                        // // b.1.2. bargain over consumption
                        int idx_power_w = index::index4(t,iL,iA,iAp,par->T,par->num_love,par->num_A,par->num_A);
                        int iPw = sol->initial_power_idx[idx_power_w];

                        int idx_power_m = index::index4(t,iL,iAp,iA,par->T,par->num_love,par->num_A,par->num_A);
                        int iPm = sol->initial_power_idx[idx_power_m];
                       
                        // b.1.3 Value conditional on meeting partner
                        //Value for woman
                        if (iPw>=0){
                            double Aw_tot = par->grid_Aw[iA] + par->grid_Am[iAp]; 
                            int idx_interp = index::couple(t,iPw,iL,0,par);
                            val_w = tools::interp_1d(par->grid_A,par->num_A,&sol->Vw_single_to_couple[idx_interp],Aw_tot);
                        } else {
                            val_w = Vw_single_to_single;
                        }

                        // Value for man
                        if (iPm>=0){
                            double Am_tot = par->grid_Am[iA] + par->grid_Aw[iAp]; 
                            int idx_interp = index::couple(t,iPm,iL,0,par);
                            val_m = tools::interp_1d(par->grid_A,par->num_A,&sol->Vm_single_to_couple[idx_interp],Am_tot);
                        } else {
                            val_m = Vm_single_to_single;
                        }

                        // expected value conditional on meeting a partner
                        Ev_cond_w += prob_w*val_w;
                        Ev_cond_m += prob_m*val_m;

                    } // iAp
                } // love 

                // b.2. expected value of starting single
                double p_meet = par->prob_repartner[t]; 
                double Ev_w = p_meet*Ev_cond_w + (1.0-p_meet)*Vw_single_to_single;
                double Ev_m = p_meet*Ev_cond_m + (1.0-p_meet)*Vm_single_to_single;

                sol->EVw_start_as_single[idx_single] = Ev_w;
                sol->EVm_start_as_single[idx_single] = Ev_m;

            } // iA
        } // pragma

        if (par->do_egm){
            calc_marginal_value_single(t, woman, sol, par);
            calc_marginal_value_single(t, man, sol, par);
        }
    }
}