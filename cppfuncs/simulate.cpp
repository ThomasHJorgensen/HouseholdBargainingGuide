
// functions for solving model for singles.
#ifndef MAIN
#define SIMULATE
#include "myheader.cpp"
#endif

namespace sim {
    void model(sim_struct *sim, sol_struct *sol, par_struct *par){
    
        // pre-compute intra-temporal optimalallocation
        #pragma omp parallel num_threads(par->threads)
        {
            #pragma omp for
            for (int i=0; i<par->simN; i++){
                for (int t=0; t < par->simT; t++){
                    int it = index::index2(i,t,par->simN,par->simT);
                    int it_1 = index::index2(i,t-1,par->simN,par->simT);

                    // state variables
                    double A_lag = sim->init_A[i];
                    double Aw_lag = sim->init_Aw[i];
                    double Am_lag = sim->init_Am[i];
                    bool couple_lag = sim->init_couple[i];
                    int power_idx_lag = sim->init_power_idx[i];
                    double love_lag = sim->init_love[i];

                    if (t>0){
                        A_lag = sim->A[it_1];
                        Aw_lag = sim->Aw[it_1];
                        Am_lag = sim->Am[it_1];
                        couple_lag = sim->couple[it_1];
                        power_idx_lag = sim->power_idx[it_1];
                        love_lag = sim->love[it_1];
                    }

                    // first check if they want to remain together and what the bargaining power will be if they do.
                    int power_idx;
                    int idx_sol = index::index4(t,power_idx_lag,0,0,par->T,par->num_power,par->num_love,par->num_A); 
                    if (couple_lag) {
                        // love-shock
                        sim->love[it] = love_lag + par->sigma_love*sim->draw_love[it];

                        // use the power index:
                        power_idx = (int) round(tools::interp_2d_int(par->grid_love,par->grid_A,par->num_love,par->num_A, &sol->power_idx[idx_sol] ,sim->love[it],A_lag));

                        if (power_idx < 0.0) { // divorce is coded as -1
                            sim->couple[it] = false;

                        } else {
                            sim->couple[it] = true;
                        }

                    } else { // remain single
                        sim->couple[it] = false;
                    }

                    // update behavior
                    if (sim->couple[it]){
                        
                        // optimal consumption allocation if couple
                        double C_tot = tools::interp_2d(par->grid_love,par->grid_A,par->num_love,par->num_A ,&sol->C_tot_couple[idx_sol],sim->love[it],A_lag);

                        double C_pub = 0.0;
                        couple::intraperiod_allocation(&sim->Cw_priv[it], &sim->Cm_priv[it], &C_pub,  C_tot,power_idx_lag,sol,par);
                        sim->Cw_pub[it] = C_pub;
                        sim->Cm_pub[it] = C_pub;

                        // update end-of-period states
                        double M_resources = par->R*A_lag + par->inc_w + par->inc_m;
                        sim->A[it] = M_resources - sim->Cw_priv[it] - sim->Cm_priv[it] - C_pub;

                        // in case of divorce
                        sim->Aw[it] = par->div_A_share * sim->A[it];
                        sim->Am[it] = (1.0-par->div_A_share) * sim->A[it];

                        sim->power_idx[it] = power_idx;
                        sim->power[it] = par->grid_power[sim->power_idx[it]];

                    } else { // single

                        // optimal consumption allocations
                        double *grid_Aw = &par->grid_A_single[index::index2(woman-1,0,2,par->num_A)];
                        double *grid_Am = &par->grid_A_single[index::index2(man-1,0,2,par->num_A)];
                        double Cw_tot = tools::interp_1d(grid_Aw,par->num_A,&sol->Cw_tot_single[t],Aw_lag);
                        double Cm_tot = tools::interp_1d(grid_Am,par->num_A,&sol->Cm_tot_single[t],Am_lag);

                        sim->Cw_priv[it] = single::cons_priv_single(Cw_tot,woman,par);
                        sim->Cw_pub[it] = Cw_tot - sim->Cw_priv[it];
                        
                        sim->Cm_priv[it] = single::cons_priv_single(Cm_tot,man,par);
                        sim->Cm_pub[it] = Cm_tot - sim->Cm_priv[it];

                        // update end-of-period states
                        double Mw = par->R*Aw_lag + par->inc_w;
                        double Mm = par->R*Am_lag + par->inc_m;
                        sim->Aw[it] = Mw - sim->Cw_priv[it] - sim->Cw_pub[it];
                        sim->Am[it] = Mm - sim->Cm_priv[it] - sim->Cm_pub[it];

                        sim->power_idx[it] = -1;

                        // left as nans by not updating them:
                        // sim->power[it] = nan
                        // sim->love[it] = nan
                        // sim->A[it] = nan
                    }

                } // t
            } // i

        } // pragma

    } // simulate
}