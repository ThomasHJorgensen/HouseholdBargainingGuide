
// functions for solving model for singles.
#ifndef MAIN
#define SIMULATE
#include "myheader.cpp"
#endif

namespace sim {

    double update_power(int t, double power_lag, double love,double A_lag,double Aw_lag,double Am_lag,sim_struct* sim, sol_struct* sol, par_struct* par){
        
        // a. value of remaining a couple at current power
        double power = 1000.0; // nonsense value
        int idx_sol = index::couple(t,0,0,0,par); 
        double Vw_couple_to_couple=0.0;
        double Vm_couple_to_couple=0.0;
        tools::interp_3d_2out(par->grid_power,par->grid_love,par->grid_A, par->num_power,par->num_love,par->num_A, &sol->Vw_couple_to_couple[idx_sol],&sol->Vm_couple_to_couple[idx_sol], power_lag,love,A_lag, &Vw_couple_to_couple, &Vm_couple_to_couple);

        // b. value of transitioning into singlehood
        int idx_single = index::single(t,0,par);
        double Vw_couple_to_single = tools::interp_1d(par->grid_Aw,par->num_A,&sol->Vw_couple_to_single[idx_single],Aw_lag);
        double Vm_couple_to_single = tools::interp_1d(par->grid_Am,par->num_A,&sol->Vm_couple_to_single[idx_single],Am_lag);
        
        // c. check participation constraints
        if ((Vw_couple_to_couple>=Vw_couple_to_single) & (Vm_couple_to_couple>=Vm_couple_to_single)){
            power = power_lag;

        } else if ((Vw_couple_to_couple<Vw_couple_to_single) & (Vm_couple_to_couple<Vm_couple_to_single)){
            power = -1.0;

        } else {
          
            // i. determine which partner is unsatisfied
            double* V_power_vec = new double[par->num_power];
            double* V_couple_to_couple = nullptr;
            double* V_couple_to_couple_partner = nullptr;
            double V_couple_to_single = 0.0;
            double V_couple_to_single_partner = 0.0;
            double* grid_power = nullptr;
            bool flip = false;
            if ((Vw_couple_to_couple<Vw_couple_to_single)){ // woman wants to leave
                V_couple_to_single = Vw_couple_to_single;
                V_couple_to_single_partner = Vm_couple_to_single;

                V_couple_to_couple = sol->Vw_couple_to_couple;
                V_couple_to_couple_partner = sol->Vm_couple_to_couple;  

                flip = false;
                grid_power = par->grid_power;                

            } else { // man wants to leave
                V_couple_to_single = Vm_couple_to_single;
                V_couple_to_single_partner = Vw_couple_to_single;

                V_couple_to_couple = sol->Vm_couple_to_couple;
                V_couple_to_couple_partner = sol->Vw_couple_to_couple;

                flip = true;
                grid_power = par->grid_power_flip;
            }

            // ii. find indifference point of unsatisfied partner:
            int j_love = tools::binary_search(0,par->num_love,par->grid_love,love); 
            int j_A = tools::binary_search(0,par->num_A,par->grid_A,A_lag); 
            for (int iP=0; iP<par->num_power; iP++){ 
                int idx = 0;
                if(flip){
                    idx = index::index4(t,par->num_power-1 - iP,0,0,par->T,par->num_power,par->num_love,par->num_A); // flipped for men
                } else {
                    idx = index::index4(t,iP,0,0,par->T,par->num_power,par->num_love,par->num_A); 
                }
                V_power_vec[iP] = tools::_interp_2d(par->grid_love,par->grid_A,par->num_love,par->num_A,&V_couple_to_couple[idx],love,A_lag,j_love,j_A);
            }
            
            // iii. interpolate the power based on the value of single to find indifference-point. (flip the axis)
            power = tools::interp_1d(V_power_vec, par->num_power, grid_power, V_couple_to_single);
            delete[] V_power_vec;
            V_power_vec = nullptr;

            if((power<0.0)|(power>1.0)){ // divorce
                power = -1.0;
            }
            else{
                // iv. find marital surplus of partner at this new power allocation
                int j_power = tools::binary_search(0,par->num_power,par->grid_power,power);
                double V_power_partner = tools::_interp_3d(par->grid_power,par->grid_love,par->grid_A, par->num_power,par->num_love,par->num_A, &V_couple_to_couple_partner[idx_sol], power,love,A_lag,j_power,j_love,j_A);
                double S_partner = couple::calc_marital_surplus(V_power_partner,V_couple_to_single_partner,par);
                
                // v. check if partner is happy. If not divorce
                if(S_partner<0.0){
                    power = -1.0; 
                }
            } 
        }

        return power;
    } // update_power


    double draw_partner_assets_cont(double A, int gender, int i, int t, sim_struct *sim, par_struct *par){
        // unpack
        double* cdf_partner_A = par->cdf_partner_Aw;
        double* grid_A = par->grid_Aw;
        double* uniform_partner_A = sim->draw_uniform_partner_Aw;
        if (gender == man){
            cdf_partner_A = par->cdf_partner_Am;
            grid_A = par->grid_Am;
            uniform_partner_A = sim->draw_uniform_partner_Am;
        }

        double* cdf_Ap_cond = new double[par->num_A];
        int index_iA = tools::binary_search(0,par->num_A,grid_A,A);

        // a. find cdf of partner assets
        for (int iAp=0; iAp<par->num_A; iAp++){
            cdf_Ap_cond[iAp] = tools::interp_1d_index_delta(grid_A,par->num_A,cdf_partner_A,A, index_iA,par->num_A, iAp,1,0);
        }

        // b. find inverted cdf of random uniform draw
        int index_sim = index::index2(i,t,par->simN,par->simT);
        double random = uniform_partner_A[index_sim];
        double A_sim = tools::interp_1d(cdf_Ap_cond,par->num_A,grid_A,random);

        delete[] cdf_Ap_cond;
        cdf_Ap_cond = nullptr;

        if (A_sim<0.0){ // WATCH OUT FOR EXTRAPOLATION OR FLAT CDF'S!!!
            A_sim = 0.0;
        }

        return A_sim;
    }

    double draw_partner_assets(double A, int gender, int i, int t, sim_struct *sim, par_struct *par){
        // unpack
        double* cdf_partner_A = par->cdf_partner_Aw;
        double* grid_A = par->grid_Aw;
        double* uniform_partner_A = sim->draw_uniform_partner_Aw;
        if (gender == man){
            cdf_partner_A = par->cdf_partner_Am;
            grid_A = par->grid_Am;
            uniform_partner_A = sim->draw_uniform_partner_Am;
        }

        // a. random uniform number
        int index_sim = index::index2(i,t,par->simN,par->simT);
        double random = uniform_partner_A[index_sim];

        // b. find first index in asset cdf above uniform draw.
        int index_iA = tools::binary_search(0,par->num_A,grid_A,A);
        for (int iAp=0; iAp<par->num_A; iAp++){
            double cdf_Ap_cond = tools::interp_1d_index_delta(grid_A,par->num_A,cdf_partner_A,A, index_iA,par->num_A, iAp,1,0);
            if(cdf_Ap_cond >= random){
                return grid_A[iAp];
            }
        }

        // c. return asset value
        return grid_A[par->num_A-1];

    }

    void model(sim_struct *sim, sol_struct *sol, par_struct *par){
    
        // pre-compute intra-temporal optimalallocation
        #pragma omp parallel num_threads(par->threads)
        {
            #pragma omp for
            for (int i=0; i<par->simN; i++){
                for (int t=0; t < par->simT; t++){
                    int it = index::index2(i,t,par->simN,par->simT);

                    // state variables
                    double A_lag = 0.0;
                    double Aw_lag = 0.0;
                    double Am_lag = 0.0;
                    bool   couple_lag = false;
                    double power_lag = 0.0;
                    double love = 0.0;
                    if (t==0){
                        A_lag = sim->init_A[i];
                        Aw_lag = sim->init_Aw[i];
                        Am_lag = sim->init_Am[i];
                        couple_lag = sim->init_couple[i];
                        power_lag = par->grid_power[sim->init_power_idx[i]];
                        love = sim->init_love[i];
                        sim->love[it] = love;
                    } else {
                        int it_1 = index::index2(i,t-1,par->simN,par->simT);
                        A_lag = sim->A[it_1];
                        Aw_lag = sim->Aw[it_1];
                        Am_lag = sim->Am[it_1];
                        couple_lag = sim->couple[it_1];
                        power_lag = sim->power[it_1];
                        love = sim->love[it];
                    } 
                    
                    // i) Find transitions in couple/single status and calculate power 
                    double power = 1000.0; // nonsense value
                    if (couple_lag) { // if start as couple

                        power = update_power(t,power_lag,love,A_lag,Aw_lag,Am_lag,sim,sol,par);
        
                        if (power < 0.0) { // divorce is coded as -1
                            sim->couple[it] = false;
                        } else {
                            sim->couple[it] = true;
                        }

                    } else { // if start as single - follow woman only
                        bool meet = (sim->draw_meet[it] < par->prob_repartner[t]);
                        if (meet){ // if meet a potential partner
                            double Ap = draw_partner_assets(Aw_lag, woman, i,t, sim, par);
                            love = sim->draw_repartner_love[it]; // note: love draws on grid.

                            power = single::calc_initial_bargaining_weight(t, love, Aw_lag, Ap, sol, par);

                            if (0.0 <= power) { // if meet and agree to couple
                                sim->couple[it] = true;

                                // set beginning-of-period couple states
                                A_lag = Aw_lag + Ap;
                                sim->love[it] = love;
                            } else { // if meet but do not agree to couple
                                power = -1.0;
                                sim->couple[it] = false;
                            }
                            
                        } else { // if do not meet
                            power = -1.0;
                            sim->couple[it] = false;
                        }
                    }

                    // ii) Find choices and update states
                    if (sim->couple[it]){
                        
                        // total consumption
                        int idx_sol = index::index4(t,0,0,0,par->T,par->num_power,par->num_love,par->num_A);
                        double C_tot = tools::interp_3d(par->grid_power,par->grid_love,par->grid_A,par->num_power,par->num_love,par->num_A ,&sol->C_tot_couple_to_couple[idx_sol],power,love,A_lag);
                        double M_resources = couple::resources(A_lag,par); // enforce ressource constraint (may be slightly broken due to approximation error)
                        if (C_tot > M_resources){ 
                            C_tot = M_resources;
                        }
                        sim->C_tot[it] = C_tot;

                        // consumpton allocation
                        double C_pub = 0.0; // placeholder for public consumption
                        precompute::intraperiod_allocation_sim(&sim->Cw_priv[it], &sim->Cm_priv[it], &C_pub,  C_tot,power,sol,par); 
                        sim->Cw_pub[it] = C_pub;
                        sim->Cm_pub[it] = C_pub;

                        // update end-of-period states
                        sim->A[it] = M_resources - sim->Cw_priv[it] - sim->Cm_priv[it] - C_pub;
                        sim->Aw[it] = par->div_A_share * sim->A[it];
                        sim->Am[it] = (1.0-par->div_A_share) * sim->A[it];
                        sim->power[it] = power;
                        if(t<par->simT-1){
                            int it1 = index::index2(i,t+1,par->simN,par->simT);
                            sim->love[it1] = love + par->sigma_love*sim->draw_love[it1];
                        }


                    } else { // single

                        // pick relevant solution for single, depending on whether just became single
                        int idx_sol_single = index::index2(t,0,par->T,par->num_A);
                        double *sol_single_w = &sol->Cw_tot_couple_to_single[idx_sol_single];
                        double *sol_single_m = &sol->Cm_tot_couple_to_single[idx_sol_single];
                        if (power_lag<0.0){
                            sol_single_w = &sol->Cw_tot_single_to_single[idx_sol_single];
                            sol_single_m = &sol->Cm_tot_single_to_single[idx_sol_single];
                        } 

                        // total consumption
                        double Cw_tot = tools::interp_1d(par->grid_Aw,par->num_A,sol_single_w,Aw_lag);
                        double Cm_tot = tools::interp_1d(par->grid_Am,par->num_A,sol_single_m,Am_lag);
                        double Mw = single::resources(Aw_lag, woman, par); // enforce ressource constraint (may be slightly broken due to approximation error)
                        double Mm = single::resources(Am_lag, man, par);
                        if (Cw_tot > Mw){
                            Cw_tot = Mw;
                        }
                        if (Cm_tot > Mm){
                            Cm_tot = Mm;
                        }
                        sim->Cm_tot[it] = Cm_tot;
                        sim->Cw_tot[it] = Cw_tot;
                        
                        // consumpton allocation
                        single::intraperiod_allocation(&sim->Cw_priv[it],&sim->Cw_pub[it],Cw_tot,woman,par);
                        single::intraperiod_allocation(&sim->Cm_priv[it],&sim->Cm_pub[it],Cm_tot,man,par);

                        // update end-of-period states  
                        sim->Aw[it] = Mw - sim->Cw_priv[it] - sim->Cw_pub[it];
                        sim->Am[it] = Mm - sim->Cm_priv[it] - sim->Cm_pub[it];
                        sim->power[it] = -1.0;

                    }

                    // iii) utility of women
                    double love_now = 0.0;
                    if (sim->couple[it]){
                        love_now = sim->love[it];
                    }
                    sim->util[it] = pow(par->beta , t) * utils::util(sim->Cw_priv[it],sim->Cw_pub[it],woman,par,love_now);

                } // t
            } // i

        } // pragma

    } // simulate
}
