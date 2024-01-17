
// functions for solving model for singles.
#ifndef MAIN
#define SIMULATE
#include "myheader.cpp"
#endif

namespace sim {

    int update_power_index(int t, int power_idx_lag, double love,double A_lag,double Aw_lag,double Am_lag,sim_struct* sim, sol_struct* sol, par_struct* par){
        
        // setup
        int power_idx = -1; // initialize as divorce
        int idx_sol = index::index4(t,power_idx_lag,0,0,par->T,par->num_power,par->num_love,par->num_A); 
        double Vw_couple, Vm_couple;

        // value of transitioning into singlehood
        int idx_single = index::index2(t,0,par->T,par->num_A);
        double Vw_single = tools::interp_1d(par->grid_Aw,par->num_A,&sol->Vw_trans_single[idx_single],Aw_lag);
        double Vm_single = tools::interp_1d(par->grid_Am,par->num_A,&sol->Vm_trans_single[idx_single],Am_lag);

        tools::interp_2d_2out(par->grid_love,par->grid_A,par->num_love,par->num_A,&sol->Vw_remain_couple[idx_sol],&sol->Vm_remain_couple[idx_sol],love,A_lag, &Vw_couple, &Vm_couple);
        
        // check participation constraints
        if ((Vw_couple>=Vw_single) & (Vm_couple>=Vm_single)){
            power_idx = power_idx_lag;

        } else {

            if ((Vm_couple>=Vm_single)){ // woman want to leave

                for (int iP=power_idx_lag+1; iP<par->num_power; iP++){ // increase power of women
                    int idx = index::index4(t,iP,0,0,par->T,par->num_power,par->num_love,par->num_A); 
                    tools::interp_2d_2out(par->grid_love,par->grid_A,par->num_love,par->num_A,&sol->Vw_remain_couple[idx],&sol->Vm_remain_couple[idx],love,A_lag, &Vw_couple, &Vm_couple);
                    
                    // check participation constraint
                    double Sw = couple::calc_marital_surplus(Vw_couple,Vw_single,par);
                    double Sm = couple::calc_marital_surplus(Vm_couple,Vm_single,par);

                    // update index if a solution is found
                    if((Sw>=0.0)&(Sm>=0.0)){
                        power_idx = iP; 
                        break;
                    }
                    
                }

            } else { // man want to leave

                for (int iP=power_idx_lag-1; iP>=0; iP--){ // increase power of men
                    int idx = index::index4(t,iP,0,0,par->T,par->num_power,par->num_love,par->num_A); 
                    tools::interp_2d_2out(par->grid_love,par->grid_A,par->num_love,par->num_A,&sol->Vw_remain_couple[idx],&sol->Vm_remain_couple[idx],love,A_lag, &Vw_couple, &Vm_couple);

                    // check participation constraint
                    double Sw = couple::calc_marital_surplus(Vw_couple,Vw_single,par);
                    double Sm = couple::calc_marital_surplus(Vm_couple,Vm_single,par);

                    // update index if a solution is found
                    if((Sw>=0.0)&(Sm>=0.0)){
                        power_idx = iP; 
                        break;
                    }
                    
                }

            }
        }

        return power_idx;
    } // update_power_index


    void model(sim_struct *sim, sol_struct *sol, par_struct *par){
    
        // pre-compute intra-temporal optimalallocation
        #pragma omp parallel num_threads(par->threads)
        {
            #pragma omp for
            for (int i=0; i<par->simN; i++){
                for (int t=0; t < par->simT; t++){
                    int it = index::index2(i,t,par->simN,par->simT);
                    int it_1 = index::index2(i,t-1,par->simN,par->simT);
                    int it1 = index::index2(i,t+1,par->simN,par->simT);

                    // state variables
                    double A_lag = sim->init_A[i];
                    double Aw_lag = sim->init_Aw[i];
                    double Am_lag = sim->init_Am[i];
                    bool couple_lag = sim->init_couple[i];
                    int power_idx_lag = sim->init_power_idx[i];
                    double love = sim->init_love[i];

                    if (t>0){
                        A_lag = sim->A[it_1];
                        Aw_lag = sim->Aw[it_1];
                        Am_lag = sim->Am[it_1];
                        couple_lag = sim->couple[it_1];
                        power_idx_lag = sim->power_idx[it_1];
                        love = sim->love[it];
                    } else {
                        sim->love[it] = love;
                    }

                    
                    // first check if they want to remain together and what the bargaining power will be if they do.
                    int power_idx;
                    if (couple_lag) {

                        power_idx = update_power_index(t,power_idx_lag,love,A_lag,Aw_lag,Am_lag,sim,sol,par);

                        if (power_idx < 0) { // divorce is coded as -1
                            sim->couple[it] = false;

                        } else {
                            sim->couple[it] = true;
                        }

                    } else { // remain single
                        sim->couple[it] = false;
                    }

                    // update behavior
                    int idx_sol = index::index4(t,power_idx,0,0,par->T,par->num_power,par->num_love,par->num_A); 
                    if (sim->couple[it]){
                        
                        // optimal consumption allocation if couple
                        double C_tot = tools::interp_2d(par->grid_love,par->grid_A,par->num_love,par->num_A ,&sol->C_tot_couple[idx_sol],love,A_lag);

                        double C_pub = 0.0;
                        couple::intraperiod_allocation(&sim->Cw_priv[it], &sim->Cm_priv[it], &C_pub,  C_tot,power_idx,sol,par); 
                        sim->Cw_pub[it] = C_pub;
                        sim->Cm_pub[it] = C_pub;

                        // update end-of-period states
                        double M_resources = couple::resources(A_lag,par); 
                        sim->A[it] = M_resources - sim->Cw_priv[it] - sim->Cm_priv[it] - C_pub;
                        if(t<par->simT-1){
                            sim->love[it1] = love + par->sigma_love*sim->draw_love[it1];
                        }

                        // in case of divorce
                        sim->Aw[it] = par->div_A_share * sim->A[it];
                        sim->Am[it] = (1.0-par->div_A_share) * sim->A[it];

                        sim->power_idx[it] = power_idx;
                        sim->power[it] = par->grid_power[power_idx];

                    } else { // single

                        // pick relevant solution for single, depending on whether just became single
                        int idx_sol_single = index::index2(t,0,par->T,par->num_A);
                        double *sol_single_w = &sol->Cw_tot_trans_single[idx_sol_single];
                        double *sol_single_m = &sol->Cm_tot_trans_single[idx_sol_single];
                        if (power_idx_lag<0){
                            sol_single_w = &sol->Cw_tot_single[idx_sol_single];
                            sol_single_m = &sol->Cm_tot_single[idx_sol_single];
                        } 

                        // optimal consumption allocations
                        double Cw_tot = tools::interp_1d(par->grid_Aw,par->num_A,sol_single_w,Aw_lag);
                        double Cm_tot = tools::interp_1d(par->grid_Am,par->num_A,sol_single_m,Am_lag);

                        single::intraperiod_allocation(&sim->Cw_priv[it],&sim->Cw_pub[it],Cw_tot,woman,par);
                        single::intraperiod_allocation(&sim->Cm_priv[it],&sim->Cm_pub[it],Cm_tot,man,par);

                        // update end-of-period states
                        double Mw = single::resources(Aw_lag,woman,par); 
                        double Mm = single::resources(Am_lag,man,par); 
                        sim->Aw[it] = Mw - sim->Cw_priv[it] - sim->Cw_pub[it];
                        sim->Am[it] = Mm - sim->Cm_priv[it] - sim->Cm_pub[it];

                        sim->power_idx[it] = -1;

                        // left as nans by not updating them:
                        // sim->power[it1] = nan
                        // sim->love[it] = nan
                        // sim->A[it] = nan
                    }

                } // t
            } // i

        } // pragma

    } // simulate
}
