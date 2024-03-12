
// functions for solving model for couples.
#ifndef MAIN
#define COUPLE
#include "myheader.cpp"
#endif

namespace couple {
    
    typedef struct {
        int t;              
        int iL;             
        int iP;             
        double M;           
        double *EVw_next;    
        double *EVm_next;    

        sol_struct *sol;
        par_struct *par;

    } solver_couple_struct;

    double calc_marital_surplus(double V_couple_to_couple,double V_couple_to_single,par_struct* par){
        return V_couple_to_couple - V_couple_to_single;
    }

    double resources(double A, par_struct* par){
        return par->R*A + par->inc_w + par->inc_m;
    }

    double value_of_choice_couple_to_couple(double* Cw_priv,double* Cm_priv,double* C_pub,double* Vw,double* Vm,  double C_tot,int t,double M_resources,int iL,int iP,double* EVw_next,double* EVm_next,sol_struct *sol, par_struct *par){
        double love = par->grid_love[iL];
        double power = par->grid_power[iP];

        // current utility from consumption allocation
        precompute::intraperiod_allocation(Cw_priv, Cm_priv, C_pub , C_tot,iP,sol,par);
        Vw[0] = utils::util(*Cw_priv,*C_pub,woman,par,love); 
        Vm[0] = utils::util(*Cm_priv,*C_pub,man,par,love);

        // add continuation value
        if (t < (par->T-1)){
            double savings = M_resources - C_tot ;
            double EVw_plus = 0.0;
            double EVm_plus = 0.0;

            EVw_plus = tools::interp_1d(par->grid_A, par->num_A, EVw_next, savings);
            EVm_plus = tools::interp_1d(par->grid_A, par->num_A, EVm_next, savings);

            Vw[0] += par->beta*EVw_plus;
            Vm[0] += par->beta*EVm_plus;
        }

        
        // return
        return power*Vw[0] + (1.0-power)*Vm[0];
    }

    //////////////////
    // VFI solution //
    double objfunc_couple_to_couple(unsigned n, const double *x, double *grad, void *solver_data_in){
        // unpack
        solver_couple_struct *solver_data = (solver_couple_struct *) solver_data_in;

        double C_tot = x[0];

        int t = solver_data->t;
        int iL = solver_data->iL;
        int iP = solver_data->iP;
        double M = solver_data->M;
        double *EVw_next = solver_data->EVw_next;
        double *EVm_next = solver_data->EVm_next;

        sol_struct *sol = solver_data->sol;
        par_struct *par = solver_data->par;

        // return negative of value
        double Cw_priv,Cm_priv,C_pub,Vw,Vm;
        return - value_of_choice_couple_to_couple(&Cw_priv,&Cm_priv,&C_pub,&Vw,&Vm, C_tot,t,M,iL,iP,EVw_next,EVm_next,sol,par);
    }

    void solve_couple_to_couple(double* Cw_priv,double* Cm_priv,double* C_pub,double* Vw,double* Vm , int t,double M_resources,int iL,int iP,double* EVw_next,double* EVm_next,double starting_val,sol_struct *sol,par_struct *par){
        
        double C_tot = M_resources;
        
        if (t<(par->T-1)){ 
            // objective function
            int const dim = 1;
            double lb[dim],ub[dim],x[dim];
            
            auto opt = nlopt_create(NLOPT_LN_BOBYQA, dim); // NLOPT_LD_MMA NLOPT_LD_LBFGS NLOPT_GN_ORIG_DIRECT
            double minf=0.0;

            solver_couple_struct* solver_data = new solver_couple_struct;
            solver_data->t = t;
            solver_data->iL = iL;
            solver_data->iP = iP;
            solver_data->M = M_resources;
            solver_data->EVw_next = EVw_next;
            solver_data->EVm_next = EVm_next;

            solver_data->sol = sol;
            solver_data->par = par;
            nlopt_set_min_objective(opt, objfunc_couple_to_couple, solver_data);
                
            // bounds
            lb[0] = 1.0e-6;
            ub[0] = solver_data->M - 1.0e-6;
            nlopt_set_lower_bounds(opt, lb);
            nlopt_set_upper_bounds(opt, ub);

            // optimize
            x[0] = starting_val;
            nlopt_optimize(opt, x, &minf);
            nlopt_destroy(opt);

            C_tot = x[0];

            delete solver_data;
        }

        // implied consumption allocation (re-calculation)
        value_of_choice_couple_to_couple(Cw_priv,Cm_priv,C_pub,Vw,Vm, C_tot,t,M_resources,iL,iP,EVw_next,EVm_next,sol,par);

    }
    
    void solve_couple_to_couple_Agrid_vfi(int t, int iP, int iL, double* EVw_next, double* EVm_next,sol_struct* sol, par_struct* par){
        for (int iA=0; iA<par->num_A;iA++){
            int idx = index::couple(t,iP,iL,iA,par);

            double M_resources = resources(par->grid_A[iA],par); 

            // starting values
            double starting_val = M_resources * 0.8;
            if (iA>0){ 
                int idx_last = index::couple(t,iP,iL,iA-1,par);
                starting_val = sol->Cw_priv_couple_to_couple[idx_last] + sol->Cm_priv_couple_to_couple[idx_last] + sol->C_pub_couple_to_couple[idx_last];
            }

            // solve unconstrained problem
            solve_couple_to_couple(&sol->Cw_priv_couple_to_couple[idx], &sol->Cm_priv_couple_to_couple[idx], &sol->C_pub_couple_to_couple[idx], &sol->Vw_couple_to_couple[idx], &sol->Vm_couple_to_couple[idx]
            , t,M_resources,iL,iP,EVw_next,EVm_next,starting_val,sol,par);
            sol->C_tot_couple_to_couple[idx] = sol->Cw_priv_couple_to_couple[idx] + sol->Cm_priv_couple_to_couple[idx] + sol->C_pub_couple_to_couple[idx];

        } // wealth   
    }


    //////////////////
    // EGM solution //

    void handle_liquidity_constraint_couple_to_couple(int t, int iP, int iL, double* m_vec, double* EmargU_pd, double* C_tot, double* Cw_priv,double* Cm_priv,double* C_pub,double* Vw,double* Vm, double* EVw_next, double* EVm_next, double* V, sol_struct* sol, par_struct* par){
        // 1. Check if liquidity constraint binds
        // constraint: binding if common m is smaller than smallest m in endogenous grid 
        for (int iA=0; iA < par->num_A; iA++){
            double M_now = resources(par->grid_A[iA],par);

            if (M_now < m_vec[0]){

                // a. Set total consumption equal to resources (consume all)
                C_tot[iA] = M_now;

                // b. Calculate intra-period allocation
                double _ = value_of_choice_couple_to_couple(&Cw_priv[iA] ,&Cm_priv[iA], &C_pub[iA], &Vw[iA], &Vm[iA], C_tot[iA],t,M_now,iL,iP,EVw_next,EVm_next,sol,par);
                
                // c. Calculate value
                double power = par->grid_power[iP];
                V[iA] = power*Vw[iA] + (1-power)*Vm[iA];
            }
        }
    }

    void do_upper_envelope_couple_to_couple(int t, int iP, int iL, double* m_vec, double* c_vec, double* v_vec, double* EmargU_pd, double* C_tot, double* Cw_priv,double* Cm_priv,double* C_pub,double* Vw,double* Vm, double* EVw_next, double* EVm_next, double* V, sol_struct* sol, par_struct* par){

        // Loop through unsorted endogenous grid
        for (int iA_pd = 0; iA_pd<par->num_A_pd-1;iA_pd++){

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
            double c_slope = (c_high - c_low)/(m_high - m_low);

            // 3. Loop through common grid
            for (int iA = 0; iA<par->num_A; iA++){

                // i. Check if resources from common grid are in current interval of endogenous grid
                double M_now = resources(par->grid_A[iA],par);
                bool interp = (M_now >= m_low) && (M_now <= m_high); 
                bool extrap_above = (iA_pd == par->num_A_pd-2) && (M_now>m_vec[par->num_A_pd-1]); // extrapolate above last point in endogenous grid
                if (interp || extrap_above){

                    // ii. Interpolate consumption and value
                    double c_guess = c_low + c_slope*(M_now - m_low);
                    double a_guess = M_now - c_guess;
                    double V_guess = V_low + v_slope*(a_guess - A_low);

                    // iii. Update sol if V is higher than previous guess (upper envelope)
                    if (V_guess > V[iA]){

                        // o. Update total consumption
                        C_tot[iA] = c_guess;
                        
                        // oo. Update intra-period allocation
                        double _ = value_of_choice_couple_to_couple(&Cw_priv[iA], &Cm_priv[iA], &C_pub[iA] ,&Vw[iA] ,&Vm[iA], C_tot[iA], t, M_now,iL,iP,EVw_next,EVm_next,sol,par);
                        
                        // ooo. Update value
                        V[iA] = par->grid_power[iP]*Vw[iA] + (1-par->grid_power[iP])*Vm[iA];
                    }
                }
            }
        }
    }

    void solve_couple_to_couple_Agrid_egm(int t, int iP, int iL, double* EVw_next, double* EVm_next, double* EmargV_next,sol_struct* sol, par_struct* par){
        
        // 1. Solve terminal period with VFI
        if(t==(par->T-1)){
            solve_couple_to_couple_Agrid_vfi(t,iP,iL,EVw_next,EVm_next,sol,par);

        // 2. Solve remaining periods with EGM
        } else {
            // Solve on endogenous grid
            double Cw_priv = 0.0;
            double Cm_priv = 0.0;
            double C_pub {};
            double Vw {};
            double Vm {};
            for (int iA_pd=0; iA_pd<par->num_A_pd;iA_pd++){

                // i. Unpack
                double A_next = par->grid_A_pd[iA_pd]; // assets next period
                int idx_pd = index::couple_pd(t,iP,iL,iA_pd,par);
                int idx_interp = index::index2(iP,0,par->num_power,par->num_marg_u);

                // ii. interpolate marginal utility
                sol->EmargU_pd[idx_pd] = par->beta * tools::interp_1d(par->grid_A, par->num_A, EmargV_next, A_next);

                // iii. Get total consumption by interpolation of pre-computed inverse marginal utility (coming from Euler)
                if (strcmp(par->interp_method,"numerical")==0){
                    
                    // starting values
                    double guess_Ctot = 3.0;
                    double guess_Cw_priv = guess_Ctot/3.0;
                    double guess_Cm_priv = guess_Ctot/3.0;
                    if(iA_pd>0){
                        // last found solution
                        guess_Ctot = sol->C_tot_pd[index::couple_pd(t,iP,iL,iA_pd-1,par)];
                        guess_Cw_priv = Cw_priv;
                        guess_Cm_priv = Cm_priv;

                    } else if (t<(par->T-2)) {
                        guess_Ctot = sol->C_tot_pd[index::couple_pd(t+1,iP,iL,iA_pd,par)];
                    }
                                        
                    sol->C_tot_pd[idx_pd] = precompute::inv_marg_util_couple(sol->EmargU_pd[idx_pd],iP,par,sol,guess_Ctot,guess_Cw_priv,guess_Cm_priv); // numerical inverse

                } else {
                    if(strcmp(par->interp_method,"linear")==0){
                        sol->C_tot_pd[idx_pd] = tools::interp_1d(&par->grid_marg_u_for_inv[idx_interp],par->num_marg_u,par->grid_inv_marg_u,sol->EmargU_pd[idx_pd]);
                    }

                    if (par->interp_inverse){
                        sol->C_tot_pd[idx_pd] = 1.0/sol->C_tot_pd[idx_pd];
                    }
                }
                
                // iv. Get endogenous grid points
                sol->M_pd[idx_pd] = A_next + sol->C_tot_pd[idx_pd];

                // v. Get post-choice value (also updates the intra-period allocation)
                sol->V_couple_to_couple_pd[idx_pd] = value_of_choice_couple_to_couple(&Cw_priv, &Cm_priv, &C_pub, &Vw, &Vm, sol->C_tot_pd[idx_pd],t,sol->M_pd[idx_pd],iL,iP,EVw_next,EVm_next,sol,par);
            }

            // 3. Apply upper envelope and interpolate onto common grid
            int idx_interp_pd = index::couple_pd(t,iP,iL,0,par);
            int idx_interp = index::couple(t,iP,iL,0,par);

            handle_liquidity_constraint_couple_to_couple(t, iP, iL, &sol->M_pd[idx_interp_pd], &sol->EmargU_pd[idx_interp_pd], 
                                                         &sol->C_tot_couple_to_couple[idx_interp], &sol->Cw_priv_couple_to_couple[idx_interp], 
                                                         &sol->Cm_priv_couple_to_couple[idx_interp], &sol->C_pub_couple_to_couple[idx_interp], 
                                                         &sol->Vw_couple_to_couple[idx_interp], &sol->Vm_couple_to_couple[idx_interp], 
                                                         EVw_next, EVm_next, &sol->V_couple_to_couple[idx_interp], sol, par);
            do_upper_envelope_couple_to_couple( t, iP, iL, 
                                                &sol->M_pd[idx_interp_pd], &sol->C_tot_pd[idx_interp_pd], &sol->V_couple_to_couple_pd[idx_interp_pd], 
                                                &sol->EmargU_pd[idx_interp_pd], &sol->C_tot_couple_to_couple[idx_interp], 
                                                &sol->Cw_priv_couple_to_couple[idx_interp] ,&sol->Cm_priv_couple_to_couple[idx_interp],
                                                &sol->C_pub_couple_to_couple[idx_interp],&sol->Vw_couple_to_couple[idx_interp],
                                                &sol->Vm_couple_to_couple[idx_interp], EVw_next, EVm_next, &sol->V_couple_to_couple[idx_interp], 
                                                sol, par);

        } // period check        
    }

    void calc_expected_value_couple(int t, int iP, int iL, double* Vw, double* Vm, double* EVw, double* EVm, sol_struct* sol, par_struct* par){
        double love = par->grid_love[iL];
        double power = par->grid_power[iP];
        for (int iA=0; iA<par->num_A;iA++){
            int idx_interp = index::couple(t,iP,0,iA,par);
            int delta_love = index::couple(t,iP,1,iA,par) - index::couple(t,iP,0,iA,par);
            double Eval_w = 0;
            double Eval_m = 0;        
            for (int i_love_shock=0; i_love_shock<par->num_shock_love; i_love_shock++){
                double love_next = love + par->grid_shock_love[i_love_shock];
                double weight = par->grid_weight_love[i_love_shock];
                int idx_love = tools::binary_search(0,par->num_love,par->grid_love,love_next);

                Eval_w += weight*tools::interp_1d_index_delta(par->grid_love, par->num_love, &Vw[idx_interp], love_next, idx_love, delta_love);
                Eval_m += weight*tools::interp_1d_index_delta(par->grid_love, par->num_love, &Vm[idx_interp], love_next, idx_love, delta_love);
            }

            int idx = index::couple(t,iP,iL,iA,par);
            EVw[idx] = Eval_w;
            EVm[idx] = Eval_m;
        }
    }

    void calc_marginal_value_couple(int t, int iP, int iL, double* Vw, double* Vm, double* margV, sol_struct* sol, par_struct* par){

        // Unpack power
        double power = par->grid_power[iP];

        // approximate marginal value of marriage by finite diff
        if (par->centered_gradient){
            for (int iA=1; iA<=par->num_A-2;iA++){
                // Setup indices
                int iA_plus = iA + 1;
                int iA_minus = iA - 1;

                // Calculate finite difference
                double margVw {0};
                double margVm {0};
                double denom = 1.0/(par->grid_A[iA_plus] - par->grid_A[iA_minus]);
                margVw = Vw[iA_plus]*denom - Vw[iA_minus]*denom;
                margVm = Vm[iA_plus]*denom - Vm[iA_minus]*denom;

                // Update solution
                margV[iA] = power*margVw + (1.0-power)*margVm;

                // Extrapolate gradient in last point
                if (iA == par->num_A-2){
                    margV[iA_plus] = margV[iA];
                }
            }

            // Extrapolate gradient in end points
            int i=0;
            margV[i] = (margV[i+2] - margV[i+1]) / (par->grid_A[i+2] - par->grid_A[i+1]) * (par->grid_A[i] - par->grid_A[i+1]) + margV[i+1];
            i = par->num_A-1;
            margV[i] = (margV[i-2] - margV[i-1]) / (par->grid_A[i-2] - par->grid_A[i-1]) * (par->grid_A[i] - par->grid_A[i-1]) + margV[i-1];
        }
        else {
            for (int iA=0; iA<=par->num_A-2;iA++){
                // Setup indices
                int iA_plus = iA + 1;

                // Calculate finite difference
                double margVw {0};
                double margVm {0};
                double denom = 1.0/(par->grid_A[iA_plus] - par->grid_A[iA]);
                margVw = Vw[iA_plus]*denom - Vw[iA]*denom;
                margVm = Vm[iA_plus]*denom - Vm[iA]*denom;

                // Update solution
                margV[iA] = power*margVw + (1.0-power)*margVm;

                // Extrapolate gradient in last point
                if (iA == par->num_A-2){
                    margV[iA_plus] = margV[iA];
                }
            }
        }
    }

    void solve_couple(int t,sol_struct *sol,par_struct *par){
        
        #pragma omp parallel num_threads(par->threads)
        {   
            // 1. Setup
            /// a. lists
            int num = 5;
            double** list_start_as_couple = new double*[num]; 
            double** list_couple_to_couple = new double*[num];
            double* list_couple_to_single = new double[num];             

            // b. temporary arrays
            double* Sw = new double[par->num_power];
            double* Sm = new double[par->num_power];

            // c. index struct to pass to bargaining algorithm
            index::index_couple_struct* idx_couple = new index::index_couple_struct;

            // 2. solve for values of remaining a couple
            #pragma omp for
            for (int iP=0; iP<par->num_power; iP++){

                // Get next period continuation values
                double *EVw_next = nullptr;  
                double *EVm_next = nullptr;
                double *EmargV_next = nullptr;
                // solve
                for (int iL=0; iL<par->num_love; iL++){
                    if (t<(par->T-1)){
                        int idx_next = index::couple(t+1,iP,iL,0,par);
                        EVw_next = &sol->EVw_start_as_couple[idx_next];  
                        EVm_next = &sol->EVm_start_as_couple[idx_next];
                        EmargV_next = &sol->EmargV_start_as_couple[idx_next];
                    }

                    if (par->do_egm){
                        solve_couple_to_couple_Agrid_egm(t,iP,iL,EVw_next,EVm_next, EmargV_next,sol,par); 

                    } else {
                        solve_couple_to_couple_Agrid_vfi(t,iP,iL,EVw_next,EVm_next,sol,par); 

                    }
                } // love
            } // power

            // 3. Solve for values of starting as couple (check participation constraints)
            #pragma omp for
            for (int iL=0; iL<par->num_love; iL++){    
                for (int iA=0; iA<par->num_A;iA++){
                    // i. Get indices
                    int idx_single = index::single(t,iA,par);
                    idx_couple->t = t;
                    idx_couple->iL = iL;
                    idx_couple->iA = iA;
                    idx_couple->par = par;

                    // ii Calculate marital surplus
                    for (int iP=0; iP<par->num_power; iP++){
                        int idx_tmp = index::couple(t,iP,iL,iA,par);
                        Sw[iP] = calc_marital_surplus(sol->Vw_couple_to_couple[idx_tmp],sol->Vw_couple_to_single[idx_single],par);
                        Sm[iP] = calc_marital_surplus(sol->Vm_couple_to_couple[idx_tmp],sol->Vm_couple_to_single[idx_single],par);
                    }

                    // iii. setup relevant lists 
                    int i = 0;
                    list_start_as_couple[i] = sol->Vw_start_as_couple; i++;
                    list_start_as_couple[i] = sol->Vm_start_as_couple; i++;
                    list_start_as_couple[i] = sol->Cw_priv_start_as_couple; i++;
                    list_start_as_couple[i] = sol->Cm_priv_start_as_couple; i++;
                    list_start_as_couple[i] = sol->C_pub_start_as_couple; i++; //consider having two of these, one for each spouse
                    i = 0;
                    list_couple_to_couple[i] = sol->Vw_couple_to_couple; i++;
                    list_couple_to_couple[i] = sol->Vm_couple_to_couple; i++;
                    list_couple_to_couple[i] = sol->Cw_priv_couple_to_couple; i++;
                    list_couple_to_couple[i] = sol->Cm_priv_couple_to_couple; i++;
                    list_couple_to_couple[i] = sol->C_pub_couple_to_couple; i++; //consider having two of these, one for each spouse
                    i = 0;
                    list_couple_to_single[i] = sol->Vw_couple_to_single[idx_single]; i++;
                    list_couple_to_single[i] = sol->Vm_couple_to_single[idx_single]; i++;
                    list_couple_to_single[i] = sol->Cw_priv_couple_to_single[idx_single]; i++;
                    list_couple_to_single[i] = sol->Cm_priv_couple_to_single[idx_single]; i++;
                    list_couple_to_single[i] = sol->Cw_pub_couple_to_single[idx_single]; i++; //consider having two of these, one for each spouse

                    // iv. Update solution
                    // Update solutions in list_start_as_couple
                    bargaining::check_participation_constraints(sol->power_idx, sol->power, Sw, Sm, idx_couple, list_start_as_couple, list_couple_to_couple, list_couple_to_single, num, par);
                    
                    // update C_tot_couple (Note: C__start_as_couple is nan when they divorce)
                    for (int iP=0; iP<par->num_power; iP++){
                        int idx = index::couple(t,iP,iL,iA,par);
                        if (sol->power[idx] >= 0.0){
                            sol->C_tot_start_as_couple[idx] = sol->Cw_priv_start_as_couple[idx] + sol->Cm_priv_start_as_couple[idx] + sol->C_pub_start_as_couple[idx];
                        }
                    }

                } // wealth
            } // love


            for (int iL=0; iL<par->num_love; iL++){
                for (int iP=0; iP<par->num_power; iP++){
                    // v. Update expected value
                    calc_expected_value_couple(t, iP, iL, sol->Vw_start_as_couple, sol->Vm_start_as_couple, sol->EVw_start_as_couple, sol->EVm_start_as_couple, sol, par);

                     // vi. Update marginal value
                    if (par->do_egm){
                        int idx_interp = index::couple(t,iP,iL,0,par);
                        // calc_marginal_value_couple(t, iP, iL, &sol->Vw_start_as_couple[idx_interp], &sol->Vm_start_as_couple[idx_interp], &sol->margV_start_as_couple[idx_interp], sol, par);
                        calc_marginal_value_couple(t, iP, iL, &sol->EVw_start_as_couple[idx_interp], &sol->EVm_start_as_couple[idx_interp], &sol->EmargV_start_as_couple[idx_interp], sol, par);
                    }
                } // power in finite diff
            } // love
            
            // delete pointers
            // for (int i=0; i<num;i++){
            //     delete[] list_start_as_couple[i];
            //     delete[] list_couple_to_couple[i];
            // }

            delete[] list_start_as_couple;
            delete[] list_couple_to_couple;
            delete[] list_couple_to_single;
            delete[] Sw;
            delete[] Sm;
            list_start_as_couple = nullptr;
            list_couple_to_couple = nullptr;
            list_couple_to_single = nullptr;
            Sw = nullptr;
            Sm = nullptr;

            delete idx_couple;

        } // pragma
    }


    void solve_single_to_couple(int t, sol_struct *sol, par_struct *par){
        #pragma omp parallel num_threads(par->threads)
        {
            #pragma omp for
            for (int iP=0; iP<par->num_power; iP++){
                for (int iL=0; iL< par->num_love; iL++){
                    for (int iA=0; iA<par->num_A;iA++){
                        int idx = index::couple(t,iP,iL,iA,par);

                        sol->Vw_single_to_couple[idx] = sol->Vw_couple_to_couple[idx];
                        sol->Vm_single_to_couple[idx] = sol->Vm_couple_to_couple[idx];
                        sol->Cw_priv_single_to_couple[idx] = sol->Cw_priv_couple_to_couple[idx];
                        sol->Cm_priv_single_to_couple[idx] = sol->Cm_priv_couple_to_couple[idx];
                        sol->C_pub_single_to_couple[idx] = sol->C_pub_couple_to_couple[idx];  
                        sol->Cw_tot_single_to_couple[idx] = sol->Cw_priv_single_to_couple[idx] + sol->C_pub_single_to_couple[idx];
                        sol->Cm_tot_single_to_couple[idx] = sol->Cm_priv_single_to_couple[idx] + sol->C_pub_single_to_couple[idx];
                    } // wealth
                } // love
            } // power
        } // pragma
    }
    
}
