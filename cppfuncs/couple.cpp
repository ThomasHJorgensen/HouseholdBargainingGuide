
// functions for solving model for singles.
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

    typedef struct{
    int t;
    int iL;
    int iA;
    par_struct *par; 

    int idx(int iP){
            return index::index4(t,iP,iL,iA , par->T,par->num_power,par->num_love,par->num_A); 
    }
    
    } index_couple_struct; 

    double calc_marital_surplus(double V_remain_couple,double V_trans_single,par_struct* par){
        return V_remain_couple - V_trans_single;
    }

    void intraperiod_allocation(double* Cw_priv, double* Cm_priv, double* C_pub , double C_tot,int iP,sol_struct *sol,par_struct *par){
        // interpolate pre-computed solution 
        int idx = index::index2(iP,0,par->num_power,par->num_Ctot);
        int j1 = tools::binary_search(0,par->num_Ctot,par->grid_Ctot,C_tot);

        Cw_priv[0] = tools::interp_1d_index(par->grid_Ctot,par->num_Ctot,&sol->pre_Ctot_Cw_priv[idx],C_tot,j1);
        Cm_priv[0] = tools::interp_1d_index(par->grid_Ctot,par->num_Ctot,&sol->pre_Ctot_Cm_priv[idx],C_tot,j1);
        // C_pub[0] = tools::interp_1d_index(par->grid_Ctot,par->num_Ctot,&sol->pre_Ctot_C_pub[idx],C_tot,j1);
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

            double M_resources = resources(par->grid_A[iA],par); //par->R*par->grid_A[iA] + par->inc_w + par->inc_m;

            // starting values
            double starting_val = M_resources * 0.8;
            if (iA>0){
                starting_val = sol->Cw_priv_remain_couple[idx_last] + sol->Cm_priv_remain_couple[idx_last] + sol->C_pub_remain_couple[idx_last];
            }

            // solve unconstrained problem
            solve_remain_couple(&sol->Cw_priv_remain_couple[idx], &sol->Cm_priv_remain_couple[idx], &sol->C_pub_remain_couple[idx], &sol->Vw_remain_couple[idx], &sol->Vm_remain_couple[idx]
            , t,M_resources,iL,iP,Vw_next,Vm_next,starting_val,sol,par);

        } // wealth   
    }

    //////////////////
    // EGM solution //
    // TODO: credit contstraint/m-grid
    void solve_remain_Agrid_egm(int t, int iP, int iL, double* Vw_next, double* Vm_next, double* marg_V_next,sol_struct* sol, par_struct* par){
        
        if(t==(par->T-1)){
            // call vfi in last period since identical: consume everything
            solve_remain_Agrid_vfi(t,iP,iL,Vw_next,Vm_next,sol,par);
            
            // store marginal utility by interpolating pre-computed values
            for (int iA=0; iA<par->num_A;iA++){
                int idx = index::index4(t,iP,iL,iA,par->T,par->num_power,par->num_love,par->num_A);
                int idx_interp = index::index2(iP,0,par->num_power,par->num_Ctot);

                double C_tot = sol->Cw_priv_remain_couple[idx] + sol->Cm_priv_remain_couple[idx] + sol->C_pub_remain_couple[idx];
                sol->marg_V_remain_couple[idx] = tools::interp_1d(par->grid_Ctot,par->num_Ctot,&par->grid_marg_u[idx_interp],C_tot);
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
                // TODO: check if this makes sense with par.R. this is when everything is in relation to M. I do not think it should now..
                double EmargU = par->beta*par->R*sol->EmargU_pd[idx_pd];
                int idx_interp = index::index2(iP,0,par->num_power,par->num_Ctot);
                
                sol->C_tot_pd[idx_pd] = tools::interp_1d(&par->grid_marg_u_for_inv[idx_interp],par->num_Ctot,par->grid_inv_marg_u,EmargU);

                // endogenous grid
                sol->M_pd[idx_pd] = A_next + sol->C_tot_pd[idx_pd];

            }

            // interpolate onto common beginning-of-period asset grid 
            // TODO: upper envelope! Not really working either
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

                // value of choice
                double _ = value_of_choice_couple(&sol->Cw_priv_remain_couple[idx] ,&sol->Cm_priv_remain_couple[idx],&sol->C_pub_remain_couple[idx],&sol->Vw_remain_couple[idx],&sol->Vm_remain_couple[idx],C_tot,t,M_now,iL,iP,Vw_next,Vm_next,sol,par);

                // marginal utility TODO: should this then not be without discounting and everything!? At least not the expected value.. MAybe that is where it went wrong..
                sol->marg_V_remain_couple[idx] = par->beta*par->R*tools::interp_1d(M_grid,par->num_A_pd,EmargU_grid,M_now);
                // sol->marg_V_remain_couple[idx] = tools::interp_1d(M_grid,par->num_A_pd,EmargU_grid,M_now);

                // int idx_interp_marg = index::index2(iP,0,par->num_power,par->num_Ctot);
                // sol->marg_V_remain_couple[idx] = tools::interp_1d(par->grid_Ctot,par->num_Ctot,&par->grid_marg_u[idx_interp_marg],C_tot);

            }

        } // period check
        
        
    }


    void update_bargaining(int *power_idx,double* power,double* Sw,double* Sm,index_couple_struct *idx_couple,double** list_start_as_couple_w,double** list_start_as_couple_m,double** list_remain_couple_w,double** list_remain_couple_m,double* list_trans_to_single_w,double* list_trans_to_single_m,int num,par_struct* par){
        
        // check the participation constraints. Array
        double min_Sw =tools::minf(Sw,par->num_power);
        double min_Sm =tools::minf(Sm,par->num_power);
        double max_Sw =tools::maxf(Sw,par->num_power);
        double max_Sm =tools::maxf(Sm,par->num_power);

        if ((min_Sw >= 0.0) & (min_Sm >= 0.0)) { // all values are consistent with marriage
            for (int iP=0; iP<par->num_power; iP++){

                // overwrite output for couple
                int idx = idx_couple->idx(iP);
                for (int i=0; i< num; i++){
                    list_start_as_couple_w[i][idx] = list_remain_couple_w[i][iP];
                    list_start_as_couple_m[i][idx] = list_remain_couple_m[i][iP];
                }
                power_idx[idx] = iP;
                power[idx] = par->grid_power[iP];
            }

        } else if ((max_Sw < 0.0) | (max_Sm < 0.0)){ // no value is consistent with marriage
            for (int iP=0; iP<par->num_power; iP++){

                // overwrite output for couple
                int idx = idx_couple->idx(iP);
                for (int i=0; i< num; i++){
                    list_start_as_couple_w[i][idx] = list_trans_to_single_w[i];
                    list_start_as_couple_m[i][idx] = list_trans_to_single_m[i];
                }
                power_idx[idx] = -1.0;
                power[idx] = -1;
            }

        } else { 

            // a. find lowest (highest) value with positive surplus for women (men)
            int Low_w = 1;      // in case there is no crossing, this will be the correct value
            int Low_m = par->num_power-1-1; // in case there is no crossing, this will be the correct value
            for (int iP=0; iP<par->num_power-1; iP++){ 
                if ((Sw[iP]<0) & (Sw[iP+1]>=0)){
                    Low_w = iP+1;
                }
                    
                if ((Sm[iP]>=0) & (Sm[iP+1]<0)){
                    Low_m = iP;
                }
            }

            // b. interpolate the surplus of each member at indifference points
            // women indifference
            int id = Low_w-1;
            double denom = (par->grid_power[id+1] - par->grid_power[id]);
            double ratio_w = (Sw[id+1] - Sw[id])/denom;
            double ratio_m = (Sm[id+1] - Sm[id])/denom;
            double power_at_zero_w = par->grid_power[id] - Sw[id]/ratio_w;
            double Sm_at_zero_w = Sm[id] + ratio_m*( power_at_zero_w - par->grid_power[id] );

            // men indifference
            id = Low_m;
            denom = (par->grid_power[id+1] - par->grid_power[id]);
            ratio_w = (Sw[id+1] - Sw[id])/denom;
            ratio_m = (Sm[id+1] - Sm[id])/denom;
            double power_at_zero_m = par->grid_power[id] - Sm[id]/ratio_m;
            double Sw_at_zero_m = Sw[id] + ratio_w*( power_at_zero_m - par->grid_power[id] );

            // c. update the outcomes
            for (int iP=0; iP<par->num_power; iP++){

                // index to store solution for couple 
                int idx = idx_couple->idx(iP);

                // i. woman wants to leave
                if (iP<Low_w){ 

                    // interpolate men's surplus
                    if (Sm_at_zero_w > 0){ // man happy to shift some bargaining power
                        for (int i=0; i< num; i++){
                            if (iP==0){
                                list_start_as_couple_w[i][idx] = tools::interp_1d_index(par->grid_power,par->num_power,list_remain_couple_w[i],power_at_zero_w,Low_w-1); 
                                list_start_as_couple_m[i][idx] = tools::interp_1d_index(par->grid_power,par->num_power,list_remain_couple_m[i],power_at_zero_w,Low_w-1); 
                            } else {
                                list_start_as_couple_w[i][idx] = list_start_as_couple_w[i][idx_couple->idx(0)]; // re-use that the interpolated values are identical
                                list_start_as_couple_m[i][idx] = list_start_as_couple_m[i][idx_couple->idx(0)]; // re-use that the interpolated values are identical
                            }
                        }
                        
                        power_idx[idx] = Low_w;
                        power[idx] = power_at_zero_w;

                    } else { // divorce

                        for (int i=0; i< num; i++){
                            list_start_as_couple_w[i][idx] = list_trans_to_single_w[i];
                            list_start_as_couple_m[i][idx] = list_trans_to_single_m[i];
                        }
                        power_idx[idx] = -1;
                        power[idx] = -1.0;
                    }
                
                } 

                // ii. man wants to leave
                else if (iP>Low_m){  

                    if (Sw_at_zero_m > 0){ // woman happy to shift some bargaining power
                        
                        for (int i=0; i< num; i++){
                            if (iP==(Low_m+1)){
                                list_start_as_couple_w[i][idx] = tools::interp_1d_index(par->grid_power,par->num_power,list_remain_couple_w[i],power_at_zero_m,Low_m); 
                                list_start_as_couple_m[i][idx] = tools::interp_1d_index(par->grid_power,par->num_power,list_remain_couple_m[i],power_at_zero_m,Low_m); 
                            } else {
                                list_start_as_couple_w[i][idx] = list_start_as_couple_w[i][idx_couple->idx(Low_m+1)]; // re-use that the interpolated values are identical
                                list_start_as_couple_m[i][idx] = list_start_as_couple_m[i][idx_couple->idx(Low_m+1)]; // re-use that the interpolated values are identical
                            }
                        }
                        power_idx[idx] = Low_m;
                        power[idx] = power_at_zero_m;
                        
                    } else { // divorce

                        for (int i=0; i< num; i++){
                            list_start_as_couple_w[i][idx] = list_trans_to_single_w[i];
                            list_start_as_couple_m[i][idx] = list_trans_to_single_m[i];
                        }

                        power_idx[idx] = -1;
                        power[idx] = -1.0;
                    }

                } 
                
                // iii. no-one wants to leave
                else { 

                    for (int i=0; i< num; i++){
                        list_start_as_couple_w[i][idx] = list_remain_couple_w[i][iP];
                        list_start_as_couple_m[i][idx] = list_remain_couple_m[i][iP];
                    }

                    power_idx[idx] = iP;
                    power[idx] = par->grid_power[iP];
                }
            } // iP

        } // outer check
        
    }

    void solve_couple(int t,sol_struct *sol,par_struct *par){
        
        #pragma omp parallel num_threads(par->threads)
        {
            // allocate memory to store relevant objects for the participation constraint check
            int shape_tmp = par->num_power;
            double* remain_Vw = new double[shape_tmp];
            double* remain_Vm = new double[shape_tmp];
            double* remain_Cw_priv = new double[shape_tmp];
            double* remain_Cm_priv = new double[shape_tmp];
            double* remain_C_pub = new double[shape_tmp];
            double* remain_marg_V = new double[shape_tmp];

            int num = 4;
            double** list_start_as_couple_w = new double*[num]; 
            double** list_start_as_couple_m = new double*[num]; 
            double** list_remain_couple_w = new double*[num]; 
            double** list_remain_couple_m = new double*[num]; 
            double* list_trans_to_single_w = new double[num]; 
            double* list_trans_to_single_m = new double[num];             

            double* Sw = new double[par->num_power];
            double* Sm = new double[par->num_power];

            index_couple_struct* idx_couple = new index_couple_struct;

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
                        solve_remain_Agrid_egm(t,iP,iL,Vw_next,Vm_next,marg_V_next,sol,par);

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
                    int idx_single_w = index::index2(t,iA,par->T,par->num_A);
                    int idx_single_m = index::index2(t,iA,par->T,par->num_A);
                    idx_couple->t = t;
                    idx_couple->iL = iL;
                    idx_couple->iA = iA;
                    idx_couple->par = par;

                    // setup temporary arrays with the one dimension being power
                    for (int iP=0; iP<par->num_power; iP++){
                        int idx_tmp = index::index4(t,iP,iL,iA,par->T,par->num_power,par->num_love,par->num_A);

                        remain_Vw[iP] = sol->Vw_remain_couple[idx_tmp];
                        remain_Vm[iP] = sol->Vm_remain_couple[idx_tmp];
                        remain_Cw_priv[iP] = sol->Cw_priv_remain_couple[idx_tmp];
                        remain_Cm_priv[iP] = sol->Cm_priv_remain_couple[idx_tmp];
                        remain_C_pub[iP] = sol->C_pub_remain_couple[idx_tmp];
                        remain_marg_V[iP] = sol->marg_V_remain_couple[idx_tmp];

                        // marital surplus
                        Sw[iP] = calc_marital_surplus(remain_Vw[iP],sol->Vw_single[idx_single_w],par);
                        Sm[iP] = calc_marital_surplus(remain_Vm[iP],sol->Vm_single[idx_single_m],par);
                    }

                    // setup relevant lists
                    int i = 0;
                    list_start_as_couple_w[i] = sol->Vw_couple; i++;
                    list_start_as_couple_w[i] = sol->Cw_priv_couple; i++;
                    list_start_as_couple_w[i] = sol->C_pub_couple; i++;
                    list_start_as_couple_w[i] = sol->marg_V_couple; i++;
                    i = 0;
                    list_start_as_couple_m[i] = sol->Vm_couple; i++;
                    list_start_as_couple_m[i] = sol->Cm_priv_couple; i++;
                    list_start_as_couple_m[i] = sol->C_pub_couple; i++;
                    list_start_as_couple_m[i] = sol->marg_V_couple; i++;

                    i = 0;
                    list_remain_couple_w[i] = remain_Vw; i++;
                    list_remain_couple_w[i] = remain_Cw_priv; i++;
                    list_remain_couple_w[i] = remain_C_pub; i++;
                    list_remain_couple_w[i] = remain_marg_V; i++;
                    i = 0;
                    list_remain_couple_m[i] = remain_Vm; i++;
                    list_remain_couple_m[i] = remain_Cm_priv; i++;
                    list_remain_couple_m[i] = remain_C_pub; i++;
                    list_remain_couple_m[i] = remain_marg_V; i++;

                    i = 0;
                    list_trans_to_single_w[i] = sol->Vw_single[idx_single_w]; i++;
                    list_trans_to_single_w[i] = sol->Cw_priv_single[idx_single_w]; i++;
                    list_trans_to_single_w[i] = sol->Cw_pub_single[idx_single_w]; i++; 
                    list_trans_to_single_w[i] = sol->Cw_pub_single[idx_single_w]; i++; // does not matter here since marg_V calcualted below
                    i = 0;
                    list_trans_to_single_m[i] = sol->Vm_single[idx_single_m]; i++;
                    list_trans_to_single_m[i] = sol->Cm_priv_single[idx_single_m]; i++;
                    list_trans_to_single_m[i] = sol->Cm_pub_single[idx_single_m]; i++; 
                    list_trans_to_single_m[i] = sol->Cm_pub_single[idx_single_m]; i++; // does not matter here since marg_V calcualted below

                    // update solution
                    update_bargaining(sol->power_idx,sol->power,Sw,Sm,idx_couple,list_start_as_couple_w,list_start_as_couple_m,list_remain_couple_w,list_remain_couple_m,list_trans_to_single_w,list_trans_to_single_m,num, par);

                    // calculate marginal utility in case of singlehood [update after check above] if EGM is implemented for singles, these numbers are stored elsewhere
                    if(par->do_egm){
                        for (int iP=0; iP<par->num_power; iP++){
                            int idx = index::index4(t,iP,iL,iA,par->T,par->num_power,par->num_love,par->num_A);
                            if (sol->power[idx] < 0.0){ // single
                                double power = par->grid_power[iP];
                                double share = par->div_A_share;

                                double Cw = sol->Cw_priv_single[idx_single_w] + sol->Cw_pub_single[idx_single_w];
                                double Cm = sol->Cm_priv_single[idx_single_m] + sol->Cm_pub_single[idx_single_m];
                                double margUw = single::marg_util_C(Cw,woman,par);
                                double margUm = single::marg_util_C(Cm,man,par);
                                sol->marg_V_couple[idx] = power*share*margUw + (1.0-power)*(1.0-share)*margUm;
                            
                            } 
                        }
                    }
                    
                } // wealth
            } // love
            
            // delete pointers
            delete[] list_start_as_couple_w;
            delete[] list_start_as_couple_m;
            delete[] list_remain_couple_w;
            delete[] list_remain_couple_m;
            delete list_trans_to_single_w;
            delete list_trans_to_single_m;

            delete Sw;
            delete Sm;

            delete remain_Vw;
            delete remain_Vm;
            delete remain_Cw_priv;
            delete remain_Cm_priv;
            delete remain_C_pub;
            delete remain_marg_V;

        } // pragma
    }
    
}
