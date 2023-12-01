
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
        double *Vw_next;    
        double *Vm_next;    

        sol_struct *sol;
        par_struct *par;

    } solver_couple_struct;

    double calc_marital_surplus(double V_remain_couple,double V_trans_single,par_struct* par){
        return V_remain_couple - V_trans_single;
    }

    void intraperiod_allocation(double* Cw_priv, double* Cm_priv, double* C_pub , double C_tot,int iP,sol_struct *sol,par_struct *par){
        // interpolate pre-computed solution 
        int idx = index::index2(iP,0,par->num_power,par->num_Ctot); 
        int j1 = tools::binary_search(0,par->num_Ctot,par->grid_Ctot,C_tot);

        Cw_priv[0] = tools::interp_1d_index(par->grid_Ctot,par->num_Ctot,&sol->pre_Ctot_Cw_priv[idx],C_tot,j1);
        Cm_priv[0] = tools::interp_1d_index(par->grid_Ctot,par->num_Ctot,&sol->pre_Ctot_Cm_priv[idx],C_tot,j1);
        C_pub[0] = C_tot - Cw_priv[0] - Cm_priv[0];

    }

    double resources(double A, par_struct* par){
        return par->R*A + par->inc_w + par->inc_m;
    }

    EXPORT double value_of_choice_couple(double* Cw_priv,double* Cm_priv,double* C_pub,double* Vw,double* Vm,  double C_tot,int t,double M_resources,int iL,int iP,double* Vw_next,double* Vm_next,sol_struct *sol, par_struct *par){
        double love = par->grid_love[iL];
        double power = par->grid_power[iP];

        // current utility from consumption allocation
        intraperiod_allocation(Cw_priv, Cm_priv, C_pub , C_tot,iP,sol,par);
        Vw[0] = utils::util(*Cw_priv,*C_pub,woman,par,love); 
        Vm[0] = utils::util(*Cm_priv,*C_pub,man,par,love);

        // add continuation value [TODO: re-use index would speed this up since only output different!]
        if (t < (par->T-1)){
            double savings = M_resources - C_tot ;
            double EVw_plus = 0.0;
            double EVm_plus = 0.0;
            for (int iL_next = 0; iL_next < par->num_shock_love; iL_next++) {
                double love_next = love + par->grid_shock_love[iL_next];

                EVw_plus += par->grid_weight_love[iL_next] * tools::interp_2d(par->grid_love,par->grid_A ,par->num_love,par->num_A, Vw_next, love_next,savings);
                EVm_plus += par->grid_weight_love[iL_next] * tools::interp_2d(par->grid_love,par->grid_A ,par->num_love,par->num_A, Vm_next, love_next,savings);
            }
            Vw[0] += par->beta*EVw_plus;
            Vm[0] += par->beta*EVm_plus;
        }

        // return
        return power*Vw[0] + (1.0-power)*Vm[0];
    }

    //////////////////
    // VFI solution //
    double objfunc_couple(unsigned n, const double *x, double *grad, void *solver_data_in){
        // unpack
        solver_couple_struct *solver_data = (solver_couple_struct *) solver_data_in;

        double C_tot = x[0];

        int t = solver_data->t;
        int iL = solver_data->iL;
        int iP = solver_data->iP;
        double M = solver_data->M;
        double *Vw_next = solver_data->Vw_next;
        double *Vm_next = solver_data->Vm_next;

        sol_struct *sol = solver_data->sol;
        par_struct *par = solver_data->par;

        // return negative of value
        double Cw_priv,Cm_priv,C_pub,Vw,Vm;
        return - value_of_choice_couple(&Cw_priv,&Cm_priv,&C_pub,&Vw,&Vm, C_tot,t,M,iL,iP,Vw_next,Vm_next,sol,par);
    }

    void solve_remain_couple(double* Cw_priv,double* Cm_priv,double* C_pub,double* Vw,double* Vm , int t,double M_resources,int iL,int iP,double* Vw_next,double *Vm_next,double starting_val,sol_struct *sol,par_struct *par){
        
        double C_tot = M_resources;
        
        if (t<(par->T-1)){ 
            // objective function
            int dim = 1;
            double lb[1],ub[1],x[1];
            
            auto opt = nlopt_create(NLOPT_LN_BOBYQA, dim); // NLOPT_LD_MMA NLOPT_LD_LBFGS NLOPT_GN_ORIG_DIRECT
            double minf=0.0;

            solver_couple_struct* solver_data = new solver_couple_struct;
            solver_data->t = t;
            solver_data->iL = iL;
            solver_data->iP = iP;
            solver_data->M = M_resources;
            solver_data->Vw_next = Vw_next;
            solver_data->Vm_next = Vm_next;

            solver_data->sol = sol;
            solver_data->par = par;
            nlopt_set_min_objective(opt, objfunc_couple, solver_data);
                
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
        }

        // implied consumption allocation (re-calculation)
        value_of_choice_couple(Cw_priv,Cm_priv,C_pub,Vw,Vm, C_tot,t,M_resources,iL,iP,Vw_next,Vm_next,sol,par);

    }
    
    void solve_remain_Agrid_vfi(int t, int iP, int iL, double* Vw_next, double* Vm_next,sol_struct* sol, par_struct* par){
        for (int iA=0; iA<par->num_A;iA++){
            int idx = index::index4(t,iP,iL,iA,par->T,par->num_power,par->num_love,par->num_A);
            int idx_last = index::index4(t,iP,iL,iA-1,par->T,par->num_power,par->num_love,par->num_A); 

            double M_resources = resources(par->grid_A[iA],par); 

            // starting values
            double starting_val = M_resources * 0.8;
            if (iA>0){ 
                starting_val = sol->Cw_priv_remain_couple[idx_last] + sol->Cm_priv_remain_couple[idx_last] + sol->C_pub_remain_couple[idx_last];
            }

            // solve unconstrained problem
            solve_remain_couple(&sol->Cw_priv_remain_couple[idx], &sol->Cm_priv_remain_couple[idx], &sol->C_pub_remain_couple[idx], &sol->Vw_remain_couple[idx], &sol->Vm_remain_couple[idx]
            , t,M_resources,iL,iP,Vw_next,Vm_next,starting_val,sol,par);
            sol->C_tot_couple[idx] = sol->Cw_priv_remain_couple[idx] + sol->Cm_priv_remain_couple[idx] + sol->C_pub_remain_couple[idx];

        } // wealth   
    }

    //////////////////
    // EGM solution //
    void solve_remain_Agrid_egm(int t, int iP, int iL, double* Vw_next, double* Vm_next, double* marg_V_next,sol_struct* sol, par_struct* par){
        
        if(t==(par->T-1)){
            // call vfi in last period since identical: consume everything
            solve_remain_Agrid_vfi(t,iP,iL,Vw_next,Vm_next,sol,par);
            
            // store marginal utility by interpolating pre-computed values
            for (int iA=0; iA<par->num_A;iA++){
                int idx = index::index4(t,iP,iL,iA,par->T,par->num_power,par->num_love,par->num_A);
                int idx_interp = index::index2(iP,0,par->num_power,par->num_marg_u);

                double C_tot = sol->Cw_priv_remain_couple[idx] + sol->Cm_priv_remain_couple[idx] + sol->C_pub_remain_couple[idx];
                sol->marg_V_remain_couple[idx] = par->beta*par->R*tools::interp_1d(par->grid_C_for_marg_u,par->num_marg_u,&par->grid_marg_u[idx_interp],C_tot);
            }

        } else {

            for (int iA_pd=0; iA_pd<par->num_A_pd;iA_pd++){
                int idx_pd = index::index3(iP,iL,iA_pd,par->num_power,par->num_love,par->num_A_pd);

                // resources next period
                double A_next = par->grid_A_pd[iA_pd];

                // interpolate marginal utility [TODO: speed up by re-using index for A and assendig order in love.]
                sol->EmargU_pd[idx_pd] = 0.0;
                for (int iL_next = 0; iL_next < par->num_shock_love; iL_next++) {
                    double love_next = par->grid_love[iL] + par->grid_shock_love[iL_next];

                    sol->EmargU_pd[idx_pd] += par->grid_weight_love[iL_next] * tools::interp_2d(par->grid_love,par->grid_A ,par->num_love,par->num_A, marg_V_next, love_next,A_next);

                }

                // invert FOC to get consumption. Using interpolation of pre-computed inverse marginal utility
                double EmargU = sol->EmargU_pd[idx_pd]; 
                int idx_interp = index::index2(iP,0,par->num_power,par->num_marg_u);
                
                sol->C_tot_pd[idx_pd] = tools::interp_1d(&par->grid_marg_u_for_inv[idx_interp],par->num_marg_u,par->grid_inv_marg_u,EmargU);

                // endogenous grid
                sol->M_pd[idx_pd] = A_next + sol->C_tot_pd[idx_pd];

                // joint value
                double Cw_priv {};
                double Cm_priv {};
                double C_pub {};
                double Vw {};
                double Vm {};

                value_of_choice_couple(&Cw_priv, &Cm_priv, &C_pub, &Vw, &Vm, sol->C_tot_pd[idx_pd],t,sol->M_pd[idx_pd],iL,iP,Vw_next,Vm_next,sol,par);

                double power = par->grid_power[iP];
                sol->V_couple_pd[idx_pd] = power*Vw + (1-power)*Vm;

            }

            // interpolate onto common beginning-of-period asset grid 
            // First draft of an upper envelope which definitly does not work
            if (par->do_upper_env){

                // endogenous grids
                int idx_interp = index::index3(iP,iL,0,par->num_power,par->num_love,par->num_A_pd);
                double* m_vec = &sol->M_pd[idx_interp];
                double* c_vec = &sol->C_tot_pd[idx_interp];
                double* v_vec = &sol->V_couple_pd[idx_interp];

                // constraint: binding if common m is smaller than smallest m in endogenous grid (check if this holds when the endo grid bends back)
                for (int iA=0; iA < par->num_A; iA++){
                    if (par->grid_A[iA] < m_vec[0]){ //<------------ check this
                        int idx = index::index4(t,iP,iL,iA,par->T,par->num_power,par->num_love,par->num_A);

                        // consume all
                        double M_now = resources(par->grid_A[iA],par);
                        double C_tot = M_now;
                        sol->C_tot_remain_couple[idx] = C_tot;

                        // Value when constrained <-------------- should we multiply by beta*R here for consistency?
                        //sol->V_remain_couple[idx] = tools::interp_1d(par->grid_C_for_marg_u, par->num_marg_u, sol->V_couple_pd, C_tot);
                        double _ = value_of_choice_couple(&sol->Cw_priv_remain_couple[idx] ,&sol->Cm_priv_remain_couple[idx],&sol->C_pub_remain_couple[idx],&sol->Vw_remain_couple[idx],&sol->Vm_remain_couple[idx],C_tot,t,M_now,iL,iP,Vw_next,Vm_next,sol,par);
                        
                        double power = par->grid_power[iP];
                        sol->V_remain_couple[idx] = power*sol->Vw_remain_couple[idx] + (1-power)*sol->Vm_remain_couple[idx];
                        // marginal  value when constrained
                        sol->marg_V_remain_couple[idx] = par->beta*par->R*tools::interp_1d(par->grid_C_for_marg_u, par->num_marg_u, par->grid_marg_u, C_tot);
                    }
                }

                // // do upper envelope
                for (int iA_pd = 0; iA_pd<par->num_A_pd-1;iA_pd++){
                    double A_low = par->grid_A_pd[iA_pd];
                    double A_high = par->grid_A_pd[iA_pd+1];

                    if (A_low>A_high){ // AMO: why is this here? would the post decision grid ever not be monotone?
                        continue;
                    }

                    // V interval and v slope AMO: for interpolation later - maybe not  necessary?
                    double V_low = v_vec[iA_pd];
                    double V_high = v_vec[iA_pd+1];
                    double v_slope = (V_high - V_low)/(A_high - A_low);

                    // m interval and c slope
                    double m_low = m_vec[iA_pd];
                    double m_high = m_vec[iA_pd+1];
                    double c_low = c_vec[iA_pd];
                    double c_high = c_vec[iA_pd+1];
                    double c_slope = (c_high - c_low)/(m_high - m_low);

                    // loop through common grid
                    for (int iA = 0; iA<par->num_A; iA++){
                        int idx = index::index4(t,iP,iL,iA,par->T,par->num_power,par->num_love,par->num_A);

                        // current m
                        double M_now = resources(par->grid_A[iA],par); 

                        // interpolate?
                        bool interp = (M_now > m_low) && (M_now < m_high); // AMO: Assumes ascending order of m_vec - DOES THIS ALWAYS HOLD?
                        //bool extrap_above = iA == par->num_A-2 && M_now>m_vec[iA_pd-1];

                        if (interp){ //} || extrap_above){

                            // implied guess
                            double c_guess = c_low + c_slope*(M_now - m_low);
                            double a_guess = M_now - c_guess;

                            // implied post decision value function
                            double V_guess = V_low + v_slope*(a_guess - A_low);

                            // update if better
                            if (V_guess > sol->V_remain_couple[idx]){
                                double C_tot = c_guess;
                                sol->C_tot_remain_couple[idx] = C_tot;
                                
                                double _ = value_of_choice_couple(&sol->Cw_priv_remain_couple[idx] ,&sol->Cm_priv_remain_couple[idx],&sol->C_pub_remain_couple[idx],&sol->Vw_remain_couple[idx],&sol->Vm_remain_couple[idx],C_tot,t,M_now,iL,iP,Vw_next,Vm_next,sol,par);
                                
                                double power = par->grid_power[iP];
                                sol->V_remain_couple[idx] = power*sol->Vw_remain_couple[idx] + (1-power)*sol->Vm_remain_couple[idx];
                                
                                // marginal value when unconstrained
                                int idx_interp = index::index3(iP,iL,0,par->num_power,par->num_love,par->num_A_pd);
                                sol->marg_V_remain_couple[idx] = par->beta*par->R*tools::interp_1d(&sol->M_pd[idx_interp],par->num_A_pd,&sol->EmargU_pd[idx_interp],M_now);


                            }

                        }

                    }
                }
            }
            else{
                for (int iA=0; iA<par->num_A;iA++){
                    int idx = index::index4(t,iP,iL,iA,par->T,par->num_power,par->num_love,par->num_A);
                    double M_now = resources(par->grid_A[iA],par); 
                    
                    // grids
                    int idx_interp = index::index3(iP,iL,0,par->num_power,par->num_love,par->num_A_pd);
                    double *M_grid = &sol->M_pd[idx_interp];
                    double *C_grid = &sol->C_tot_pd[idx_interp];
                    double *EmargU_grid = &sol->EmargU_pd[idx_interp];

                    // total consumption
                    double C_tot = tools::interp_1d(M_grid,par->num_A_pd,C_grid,M_now);

                    // check credit constraint
                    if (C_tot > M_now){
                        C_tot = M_now;

                        // marginal value when constrained
                        sol->marg_V_remain_couple[idx] = par->beta*par->R*tools::interp_1d(par->grid_C_for_marg_u, par->num_marg_u, par->grid_marg_u, C_tot); 
                    }
                    else{
                        // marginal value when unconstrained
                        sol->marg_V_remain_couple[idx] = par->beta*par->R*tools::interp_1d(M_grid,par->num_A_pd,EmargU_grid,M_now);
                    }

                    // value of choice (allocate consumption)
                    double _ = value_of_choice_couple(&sol->Cw_priv_remain_couple[idx] ,&sol->Cm_priv_remain_couple[idx],&sol->C_pub_remain_couple[idx],&sol->Vw_remain_couple[idx],&sol->Vm_remain_couple[idx],C_tot,t,M_now,iL,iP,Vw_next,Vm_next,sol,par);

                } //common grid
            } // upper envelope

        } // period check
        
        
    }


    void solve_couple(int t,sol_struct *sol,par_struct *par){
        
        #pragma omp parallel num_threads(par->threads)
        {   int num = 6;
            double** list_start_as_couple = new double*[num]; 
            double** list_remain_couple = new double*[num];
            double* list_trans_to_single = new double[num];             

            double* Sw = new double[par->num_power];
            double* Sm = new double[par->num_power];

            index::index_couple_struct* idx_couple = new index::index_couple_struct;

            // a. solve for values of reminaing a couple
            #pragma omp for
            for (int iP=0; iP<par->num_power; iP++){

                // continuation values
                int idx_next = index::index4(t+1,iP,0,0,par->T,par->num_power,par->num_love,par->num_A);
                if (t==(par->T-1)){ // does not matter in last period-> fix at some valid index
                    idx_next = 0;
                }
                double *Vw_next = &sol->Vw_couple[idx_next];  
                double *Vm_next = &sol->Vm_couple[idx_next];
                double *marg_V_next = &sol->marg_V_couple[idx_next];

                for (int iL=0; iL<par->num_love; iL++){

                    // solve for all values in grid_A.
                    if (par->do_egm){
                        // solve_remain_Agrid_egm(t,iP,iL,Vw_next,Vm_next,marg_V_next,sol,par); 
                        solve_remain_Agrid_egm(t,iP,iL,Vw_next,Vm_next, marg_V_next,sol,par); 

                    } else {
                        solve_remain_Agrid_vfi(t,iP,iL,Vw_next,Vm_next,sol,par); 

                    }
 
                } // love
            } // power

            // b. check the participation constraints. Same loops just without iP
            #pragma omp for
            for (int iL=0; iL<par->num_love; iL++){    
                for (int iA=0; iA<par->num_A;iA++){
                    // indices
                    int idx_single = index::index2(t,iA,par->T,par->num_A);
                    idx_couple->t = t;
                    idx_couple->iL = iL;
                    idx_couple->iA = iA;
                    idx_couple->par = par;

                    // setup temporary array of marital surplus with the one dimension being power
                    for (int iP=0; iP<par->num_power; iP++){
                        int idx_tmp = index::index4(t,iP,iL,iA,par->T,par->num_power,par->num_love,par->num_A);
                        Sw[iP] = calc_marital_surplus(sol->Vw_remain_couple[idx_tmp],sol->Vw_single[idx_single],par);
                        Sm[iP] = calc_marital_surplus(sol->Vm_remain_couple[idx_tmp],sol->Vm_single[idx_single],par);
                    }

                    // setup relevant lists
                    int i = 0;
                    list_start_as_couple[i] = sol->Vw_couple; i++;
                    list_start_as_couple[i] = sol->Vm_couple; i++;
                    list_start_as_couple[i] = sol->Cw_priv_couple; i++;
                    list_start_as_couple[i] = sol->Cm_priv_couple; i++;
                    list_start_as_couple[i] = sol->C_pub_couple; i++; //consider having two of these, one for each spouse
                    list_start_as_couple[i] = sol->marg_V_couple; 
                    i = 0;
                    list_remain_couple[i] = sol->Vw_remain_couple; i++;
                    list_remain_couple[i] = sol->Vm_remain_couple; i++;
                    list_remain_couple[i] = sol->Cw_priv_remain_couple; i++;
                    list_remain_couple[i] = sol->Cm_priv_remain_couple; i++;
                    list_remain_couple[i] = sol->C_pub_remain_couple; i++; //consider having two of these, one for each spouse
                    list_remain_couple[i] = sol->marg_V_remain_couple; 
                    i = 0;
                    list_trans_to_single[i] = sol->Vw_single[idx_single]; i++;
                    list_trans_to_single[i] = sol->Vm_single[idx_single]; i++;
                    list_trans_to_single[i] = sol->Cw_priv_single[idx_single]; i++;
                    list_trans_to_single[i] = sol->Cm_priv_single[idx_single]; i++;
                    list_trans_to_single[i] = sol->Cw_pub_single[idx_single]; i++; //consider having two of these, one for each spouse
                    list_trans_to_single[i] = sol->marg_Vw_single[idx_single]; // doesn't matter because not used in EGM

                    // update solution
                    bargaining::check_participation_constraints(sol->power_idx, sol->power, Sw, Sm, idx_couple, list_start_as_couple, list_remain_couple, list_trans_to_single, num, par);

                    // calculate marginal utility in case of singlehood [update after check above] if EGM is implemented for singles, these numbers are stored elsewhere
                    if(par->do_egm){
                        for (int iP=0; iP<par->num_power; iP++){
                            int idx = index::index4(t,iP,iL,iA,par->T,par->num_power,par->num_love,par->num_A);
                            if (sol->power[idx] < 0.0){ // single
                                double power = par->grid_power[iP];
                                double share = par->div_A_share;

                                double Cw = sol->Cw_priv_single[idx_single] + sol->Cw_pub_single[idx_single];
                                double Cm = sol->Cm_priv_single[idx_single] + sol->Cm_pub_single[idx_single];
                                double margUw = utils::marg_util_C(Cw,woman,par); 
                                double margUm = utils::marg_util_C(Cm,man,par);
                                sol->marg_V_couple[idx] = power*share*margUw + (1.0-power)*(1.0-share)*margUm; 
                            
                            } 
                        }
                    }
                    
                } // wealth
            } // love
            
            // delete pointers
            delete[] list_start_as_couple;
            delete[] list_remain_couple;
            delete list_trans_to_single;

            delete Sw;
            delete Sm;

        } // pragma
    }
    
}
