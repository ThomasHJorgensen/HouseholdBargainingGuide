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


    EXPORT void handle_liquidity_constraint_single_to_single(int t, int gender, double* m_vec, double* C_tot, double* C_priv, double* C_pub, double* V, double* EV_next, sol_struct* sol, par_struct* par){
        // 1. Check if liquidity constraint binds
        // constraint: binding if common m is smaller than smallest m in endogenous grid
        double* grid_A = par->grid_Aw;
        if (gender == man) {
            grid_A = par->grid_Am;
        }
        for (int iA = 0; iA < par->num_A; iA++){
            double M_now = resources(grid_A[iA],gender,par);

            if (M_now < m_vec[0]){

                // a. Set total consumption equal to resources (consume all)
                C_tot[iA] = M_now;

                // b. Calculate intra-period allocation
                intraperiod_allocation(&C_priv[iA], &C_pub[iA], C_tot[iA], gender, par);

                // c. Calculate value
                V[iA] = value_of_choice_single_to_single(C_tot[iA], M_now, gender, EV_next, par);

            }
        }
    }

    void do_upper_envelope_single_to_single(int t, int gender, double* m_vec, double* c_vec, double* v_vec, double* C_tot, double* C_priv, double* C_pub, double* V, double* EV_next, sol_struct* sol, par_struct* par){
        // Unpack
        double* grid_A = par->grid_Aw;
        if (gender == man){
            grid_A = par->grid_Am;
        }

        // Loop through unsorted endogenous grid
        for (int iA_pd = 0; iA_pd<par->num_A_pd; iA_pd++){

            // 1. Unpack intervals
            double A_low = par->grid_A_pd[iA_pd];
            double A_high = par->grid_A_pd[iA_pd+1];
            
            double V_low = v_vec[iA_pd];
            double V_high = v_vec[iA_pd+1];

            double m_low = m_vec[iA_pd];
            double m_high = m_vec[iA_pd+1];

            double c_low = c_vec[iA_pd];
            double c_high = c_vec[iA_pd+1];

            // 2. Calculate slopes
            double v_slope = (V_high - V_low)/(A_high - A_low);
            double c_slope = (c_high - c_low)/(A_high - A_low);

            // 3. Loop through common grid
            for (int iA = 0; iA<par->num_A; iA++){

                // i. Check if resources from common grid are in current interval of endogenous grid
                double M_now = resources(grid_A[iA], gender, par);
                bool interp = ((M_now >= m_low) & (M_now <= m_high));
                bool extrap_above = ((iA_pd == par->num_A_pd-2) & (M_now > m_vec[par->num_A_pd-1])); // extrapolate above last point in endogenous grid

                if (interp | extrap_above){

                    // ii. Interpolate consumption and value
                    double c_guess = c_low + c_slope*(M_now - m_low);
                    double a_guess = M_now - c_guess;
                    double V_guess = V_low + v_slope*(a_guess - A_low);

                    // iii. Update sol if v is higher than previous guess (upper envelope)
                    if (V_guess > V[iA]){
                        // o. Update total consumption
                        C_tot[iA] = c_guess;

                        // oo. Update intra-period allocation
                        intraperiod_allocation(&C_priv[iA], &C_pub[iA], C_tot[iA], gender, par);

                        // ooo. Update value
                        V[iA] = value_of_choice_single_to_single(C_tot[iA], M_now, gender, EV_next, par);
                    }
                }
            }
        }
    }


    void EGM_single_to_single(int t, int gender, sol_struct* sol, par_struct* par){
        // 1. Setup
        /// a. unpack
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
        double* V_pd {sol->Vw_single_to_single_pd};
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
            V_pd = sol->Vm_single_to_single_pd;
        }

        /// c. Allocate memory
        double* EmargU_pd {new double[par->num_A_pd]};
        double* C_tot_pd {new double[par->num_A_pd]};
        double* M_pd {new double[par->num_A_pd]};

        // 2. EGM step
        /// setup
        int idx = index::single(t,0, par);
        int idx_next = index::single(t+1,0, par);
        int min_point_A = 0;

        for (int iA_pd=0; iA_pd<par->num_A_pd; iA_pd++){

            /// a. get next period assets
            double A_next = grid_A_pd[iA_pd];

            /// b. calculate expected marginal utility
            min_point_A = tools::binary_search(min_point_A, par->num_A, grid_A, A_next);
            EmargU_pd[iA_pd] = par->beta*tools::interp_1d_index(grid_A, par->num_A, &margV[idx_next],A_next, min_point_A);

            /// c. invert marginal utility by interpolation from pre-computed grid
            C_tot_pd[iA_pd] = utils::inv_marg_util_C(EmargU_pd[iA_pd], gender, par);

            /// d. endogenous grid over resources
            M_pd[iA_pd] = C_tot_pd[iA_pd] + A_next;

            /// e. value
            V_pd[iA_pd] = value_of_choice_single_to_single(C_tot_pd[iA_pd], M_pd[iA_pd], gender, &EV[idx_next], par);
        }

        // 3. liquidity constraint
        handle_liquidity_constraint_single_to_single(t, gender, M_pd, &C_tot[idx], &C_priv[idx], &C_pub[idx], &V[idx], &EV[idx_next], sol, par);

        // 4. upper envelope
        do_upper_envelope_single_to_single(t, gender, M_pd, C_tot_pd, V_pd, &C_tot[idx], &C_priv[idx], &C_pub[idx], &V[idx], &EV[idx_next], sol, par);


        // 4. clean up
        delete[] EmargU_pd;
        delete[] C_tot_pd;
        delete[] M_pd;
        EmargU_pd = nullptr;
        C_tot_pd = nullptr;
        M_pd = nullptr;
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
        if (par->centered_gradient){
            for (int iA=1; iA<num_A-1; iA++){
                // Setup indices
                int iA_plus = iA + 1;
                int iA_minus = iA - 1;

                double denom = 1/(grid_A[iA_plus] - grid_A[iA_minus]);

                // Calculate finite difference
                margV[iA] = V[iA_plus]*denom - V[iA_minus]* denom; 
            }
             // Extrapolate gradient in end points
            int i=0;
            margV[i] = (margV[i+2] - margV[i+1]) / (grid_A[i+2] - grid_A[i+1]) * (grid_A[i] - grid_A[i+1]) + margV[i+1];
            i = par->num_A-1;
            margV[i] = (margV[i-2] - margV[i-1]) / (grid_A[i-2] - grid_A[i-1]) * (grid_A[i] - grid_A[i-1]) + margV[i-1];
            
        } 
        else {
            for (int iA=0; iA<num_A-1; iA++){
                // Setup indices
                int iA_plus = iA + 1;

                double denom = 1/(grid_A[iA_plus] - grid_A[iA]);

                // Calculate finite difference
                margV[iA] = V[iA_plus]*denom - V[iA]* denom; 

                // Extrapolate gradient in last point
                if (iA == num_A-2){
                    margV[iA_plus] = margV[iA];
                }
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
                    delete solver_data;

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

     double repartner_surplus(double power, index::state_couple_struct* state_couple, index::state_single_struct* state_single, int gender, par_struct* par, sol_struct* sol){ //TODO: add index
        // unpack
        int t = state_single->t;
        double A = state_single->A;
        double love = state_couple->love;
        double A_tot = state_couple->A; 

        // gender specific
        double* V_single_to_single = sol->Vw_single_to_single;
        double* V_single_to_couple = sol->Vw_single_to_couple;
        double* grid_A_single = par->grid_Aw;
        if (gender == man){
            V_single_to_single = sol->Vm_single_to_single;
            V_single_to_couple = sol->Vm_single_to_couple;
            grid_A_single = par->grid_Am;
        }
        
        // Get indices
        int iA_single = state_single->iA;
        int iL_couple = state_couple->iL;
        int iA_couple = state_couple->iA;
        int iP = tools::binary_search(0, par->num_power, par->grid_power, power);

        // Get indices if not provided
        if (iL_couple == -1){
            iL_couple = tools::binary_search(0, par->num_love, par->grid_love, love);
        }
        if (iA_couple == -1){
            iA_couple = tools::binary_search(0, par->num_A, par->grid_A, A_tot);
        }
        if (iA_single == -1){
            iA_single = tools::binary_search(0, par->num_A, grid_A_single, A);
        }

        //interpolate V_single_to_single
        int idx_single = index::single(t,0,par);
        double Vsts = tools::interp_1d_index(grid_A_single, par->num_A, &V_single_to_single[idx_single], A, iA_single); 

        // interpolate couple V_single_to_couple  
        int idx_couple = index::couple(t,0,0,0,par);
        double Vstc = tools::_interp_3d(par->grid_power, par->grid_love, par->grid_A, 
                                       par->num_power, par->num_love, par->num_A, 
                                       &V_single_to_couple[idx_couple], power, love, A_tot,
                                       iP, iL_couple, iA_couple);

        // surplus
        return Vstc - Vsts;
    }

    double calc_initial_bargaining_weight(int t, double love, double Aw, double Am, sol_struct* sol, par_struct* par, int iL_couple=-1){ //TODO: add index
        // state structs
        index::state_couple_struct* state_couple = new index::state_couple_struct;
        index::state_single_struct* state_single_w = new index::state_single_struct;
        index::state_single_struct* state_single_m = new index::state_single_struct;

        // couple
        state_couple->t = t;
        state_couple->love = love;
        state_couple->A = Aw+Am;
        state_couple->iA = tools::binary_search(0, par->num_A, par->grid_A, Aw+Am);
        if (iL_couple == -1){
            iL_couple = tools::binary_search(0, par->num_love, par->grid_love, love);
        }
        state_couple->iL = iL_couple;

        // single woman
        state_single_w->t = t;
        state_single_w->A = Aw;
        state_single_w->iA = tools::binary_search(0, par->num_A, par->grid_Aw, Aw);

        // single man
        state_single_m->t = t;
        state_single_m->A = Am;
        state_single_m->iA = tools::binary_search(0, par->num_A, par->grid_Am, Am);
        // Note: We don't know whether we are on the woman or man asset grid, so we need to search both.
        // We could pass gender to calc_initial_bargaining_weight to infer which grid we are on, and avoid binary search for that gender

        //solver input
        bargaining::nash_solver_struct* nash_struct = new bargaining::nash_solver_struct;
        nash_struct->surplus_func = repartner_surplus;
        nash_struct->state_couple = state_couple;
        nash_struct->state_single_w = state_single_w;
        nash_struct->state_single_m = state_single_m;
        nash_struct->sol = sol;
        nash_struct->par = par;

        // solve
        double init_mu =  bargaining::nash_bargain(nash_struct);

        delete state_couple;
        delete state_single_w;
        delete state_single_m;
        delete nash_struct;

        return init_mu;
    }
    
    
    double expected_value_cond_meet_partner(int t, int iA, int gender, sol_struct* sol, par_struct* par){
        // unpack
        double* V_single_to_single = sol->Vw_single_to_single;
        double* V_single_to_couple = sol->Vw_single_to_couple;
        double* prob_partner_A = par->prob_partner_A_w;
        double* grid_A = par->grid_Aw;
        if (gender == man){
            V_single_to_single = sol->Vm_single_to_single;
            V_single_to_couple = sol->Vm_single_to_couple;
            prob_partner_A = par->prob_partner_A_m;
            grid_A = par->grid_Am;
        }
        // // value of remaining single
        int idx_single = index::single(t,iA,par);

        // // b.1. loop over potential partners conditional on meeting a partner
        double Ev_cond = 0.0;
        double val = 0.0;
        for(int iL=0;iL<par->num_love;iL++){
            for(int iAp=0;iAp<par->num_A;iAp++){ // partner's wealth 

                // b.1.1. probability of meeting a specific type of partner
                int idx_A = index::index2(iA,iAp,par->num_A,par->num_A);
                double prob_A = prob_partner_A[idx_A]; 
                double prob_love = par->prob_partner_love[iL]; 
                double prob = prob_A*prob_love;

                // only calculate if match has positive probability of happening
                if (prob>0.0) {
                    // Figure out gender
                    int iAw = iA;
                    int iAm = iAp;
                    if (gender==man) {
                        int iAw = iAp;
                        int iAm = iA;
                    }

                    // // b.1.2. bargain over consumption
                    double love = par->grid_love[iL];
                    double Aw = grid_A[iAw];
                    double Am = grid_A[iAm];
                    double power = calc_initial_bargaining_weight(t, love, Aw, Am, sol, par, iL);
                    
                    // b.1.3 Value conditional on meeting partner
                    if (power>=0.0){
                        double A_tot = Aw + Am;
                        int idx_interp = index::couple(t, 0, 0, 0, par);
                        val = tools::interp_3d(par->grid_power, par->grid_love, par->grid_A, 
                                       par->num_power, par->num_love, par->num_A, 
                                       &V_single_to_couple[idx_interp], power, love, A_tot); //TODO: reuse index
                    } else {
                        val = V_single_to_single[idx_single];
                    }

                    // expected value conditional on meeting a partner
                    Ev_cond += prob*val;
                } // if
            } // iAp
        } // love 
        return Ev_cond;
    }



    void expected_value_start_single(int t, sol_struct* sol,par_struct* par){
        #pragma omp parallel num_threads(par->threads)
        {// a. Loop over states
            #pragma omp for
            for (int iA=0; iA<par->num_A;iA++){
                // a.1 Value conditional on meeting partner
                double EVw_cond = expected_value_cond_meet_partner(t,iA,woman,sol,par);
                double EVm_cond = expected_value_cond_meet_partner(t,iA,man,sol,par);

                // a.2. expected value of starting single
                double p_meet = par->prob_repartner[t];
                int idx_single = index::single(t,iA,par);
                sol->EVw_start_as_single[idx_single] = p_meet*EVw_cond + (1.0-p_meet)*sol->Vw_single_to_single[idx_single];
                sol->EVm_start_as_single[idx_single] = p_meet*EVm_cond + (1.0-p_meet)*sol->Vm_single_to_single[idx_single];

                sol->EVw_cond_meet_partner[idx_single] = EVw_cond;
                sol->EVm_cond_meet_partner[idx_single] = EVm_cond;
            } // iA
        } // pragma

        if (par->do_egm){
            calc_marginal_value_single(t, woman, sol, par);
            calc_marginal_value_single(t, man, sol, par);
        }
    }

}
