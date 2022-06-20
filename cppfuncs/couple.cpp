
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

    void intraperiod_allocation(double* Cw_priv, double* Cm_priv, double* C_pub , double C_tot,int iP,sol_struct *sol,par_struct *par){
        // interpolate pre-computed solution 
        int idx = index::index2(iP,0,par->num_power,par->num_Ctot);
        int j1 = tools::binary_search(0,par->num_Ctot,par->grid_Ctot,C_tot);

        Cw_priv[0] = tools::interp_1d_index(par->grid_Ctot,par->num_Ctot,&sol->pre_Ctot_Cw_priv[idx],C_tot,j1);
        Cm_priv[0] = tools::interp_1d_index(par->grid_Ctot,par->num_Ctot,&sol->pre_Ctot_Cm_priv[idx],C_tot,j1);
        C_pub[0] = tools::interp_1d_index(par->grid_Ctot,par->num_Ctot,&sol->pre_Ctot_C_pub[idx],C_tot,j1);

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

    EXPORT void solve_uncon_couple(double* Cw_priv,double* Cm_priv,double* C_pub,double* Vw,double* Vm , int t,double M_resources,int iL,int iP,double* Vw_next,double *Vm_next,double starting_val,sol_struct *sol,par_struct *par){
        
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


    void check_participation_constraints(int *power,double* Sw,double* Sm,int idx_single,index_couple_struct *idx_couple,double** list_couple,double** list_raw,double** list_single,int num,par_struct* par){
        
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
                    list_couple[i][idx] = list_raw[i][iP];
                }
                power[idx] = iP;
            }

        } else if ((max_Sw < 0.0) | (max_Sm < 0.0)){ // no value is consistent with marriage
            for (int iP=0; iP<par->num_power; iP++){

                // overwrite output for couple
                int idx = idx_couple->idx(iP);
                for (int i=0; i< num; i++){
                    list_couple[i][idx] = list_single[i][idx_single];
                }
                power[idx] = -1;
            }

        } else { 

            // a. find lowest (highest) value with positive surplus for women (men)
            int Low_w = 0;      // in case there is no crossing, this will be the correct value
            int Low_m = par->num_power-1; // in case there is no crossing, this will be the correct value
            for (int iP=0; iP<par->num_power; iP++){
                if ((Sw[iP]<0) & (Sw[iP+1]>=0)){
                    Low_w = iP+1;
                }
                    
                if ((Sm[iP]>=0) & (Sm[iP+1]<0)){
                    Low_m = iP;
                }
            }

            // b. update the outcomes
            for (int iP=0; iP<par->num_power; iP++){

                // index to store solution for couple 
                int idx = idx_couple->idx(iP);

                // i. woman wants to leave
                if (iP<Low_w){ 
                    if (Sm[Low_w] > 0){ // man happy to shift some bargaining power

                        for (int i=0; i< num; i++){
                            list_couple[i][idx] = list_raw[i][Low_w];
                        }
                        power[idx] = Low_w;
                    } else { // divorce

                        for (int i=0; i< num; i++){
                            list_couple[i][idx] = list_single[i][idx_single];
                        }
                        power[idx] = -1;
                    }
                
                } 

                // ii. man wants to leave
                else if (iP>Low_m){  
                    if (Sw[Low_m] > 0){ // woman happy to shift some bargaining power
                        
                        for (int i=0; i< num; i++){
                            list_couple[i][idx] = list_raw[i][Low_m];
                        }
                        power[idx] = Low_m;
                        
                    } else { // divorce

                        for (int i=0; i< num; i++){
                            list_couple[i][idx] = list_single[i][idx_single];
                        }

                        power[idx] = -1;
                    }

                } 
                
                // iii. no-one wants to leave
                else { 

                    for (int i=0; i< num; i++){
                        list_couple[i][idx] = list_raw[i][iP];
                    }

                    power[idx] = iP;
                }
            } // iP

        } // outer check
        
    }


    void solve_couple(int t,sol_struct *sol,par_struct *par){
        
        #pragma omp parallel num_threads(par->threads)
        {
            // allocate memory to store relevant objects
            double* tmp_Vw = new double[par->num_power];
            double* tmp_Vm = new double[par->num_power];
            double* tmp_Cw_priv = new double[par->num_power];
            double* tmp_Cm_priv = new double[par->num_power];
            double* tmp_C_pub = new double[par->num_power];

            int num = 5;
            double** list_couple = new double*[num]; 
            double** list_single = new double*[num]; 
            double** list_raw = new double*[num]; 

            double* Sw = new double[par->num_power];
            double* Sm = new double[par->num_power];

            index_couple_struct* idx_couple = new index_couple_struct;

            #pragma omp for
            for (int iL=0; iL<par->num_love; iL++){
                double love = par->grid_love[iL];
                for (int iA=0; iA<par->num_A;iA++){
                    double M_resources = par->R*par->grid_A[iA] + par->inc_w + par->inc_m;
                    
                    double starting_val = M_resources * 0.8;
                    for (int iP=0; iP<par->num_power; iP++){
                        double power = par->grid_power[iP];

                        // continuation values
                        int idx_next = index::index4(t+1,iP,0,0,par->T,par->num_power,par->num_love,par->num_A);
                        if (t==(par->T-1)){ // does not matter in last period-> fix at some valid index
                            idx_next = 0;
                        }
                        double *Vw_next = &sol->Vw_couple[idx_next];
                        double *Vm_next = &sol->Vm_couple[idx_next];

                        // starting values
                        if (iP>0){
                            double C_tot_last = tmp_Cw_priv[iP-1] + tmp_Cm_priv[iP-1] + tmp_C_pub[iP-1];
                            starting_val = C_tot_last;
                        }

                        // solve unconstrained problem
                        solve_uncon_couple(&tmp_Cw_priv[iP], &tmp_Cm_priv[iP], &tmp_C_pub[iP], &tmp_Vw[iP], &tmp_Vm[iP] , t,M_resources,iL,iP,Vw_next,Vm_next,starting_val,sol,par);
                    
                    } // power

                    // check the participation constraints
                    // indices
                    int idx_single = index::index2(t,iA,par->T,par->num_A);
                    idx_couple->t = t;
                    idx_couple->iL = iL;
                    idx_couple->iA = iA;
                    idx_couple->par = par;

                    // setup relevant lists
                    int i = 0;
                    list_couple[i] = sol->Vw_couple; i++;
                    list_couple[i] = sol->Vm_couple; i++;
                    list_couple[i] = sol->Cw_priv_couple; i++;
                    list_couple[i] = sol->Cm_priv_couple; i++;
                    list_couple[i] = sol->C_pub_couple; i++;

                    i = 0;
                    list_raw[i] = tmp_Vw; i++;
                    list_raw[i] = tmp_Vm; i++;
                    list_raw[i] = tmp_Cw_priv; i++;
                    list_raw[i] = tmp_Cm_priv; i++;
                    list_raw[i] = tmp_C_pub; i++;

                    i = 0;
                    list_single[i] = sol->Vw_single; i++;
                    list_single[i] = sol->Vm_single; i++;
                    list_single[i] = sol->Cw_priv_single; i++;
                    list_single[i] = sol->Cm_priv_single; i++;
                    list_single[i] = sol->Cw_pub_single; i++; // does not matter
                    
                    // marital surplus
                    for (int iP=0; iP<par->num_power; iP++){
                        Sw[iP] = tmp_Vw[iP] - sol->Vw_single[idx_single];
                        Sm[iP] = tmp_Vm[iP] - sol->Vm_single[idx_single];
                    }

                    // update solution
                    check_participation_constraints(sol->power_idx,Sw,Sm,idx_single,idx_couple,list_couple,list_raw,list_single,num, par);
                    
                }
            }
            
            // delete pointers
            delete[] list_couple;
            delete[] list_single;
            delete[] list_raw;

            delete Sw;
            delete Sm;

            delete tmp_Vw;
            delete tmp_Vm;
            delete tmp_Cw_priv;
            delete tmp_Cm_priv;
            delete tmp_C_pub;

        } // pragma
    }


    
}