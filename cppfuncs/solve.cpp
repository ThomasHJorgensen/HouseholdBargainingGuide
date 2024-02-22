#define MAIN
#include "myheader.h"

// include these again here to ensure that they are automatically compiled by consav
#ifndef MAIN
#include "precompute.cpp"
#endif

/////////////
// 5. MAIN //
/////////////

EXPORT void solve(sol_struct *sol, par_struct *par){
    
    // pre-compute intra-temporal optimal allocation
    precompute::precompute(sol,par);

    // loop backwards
    for (int t = par->T-1; t >= 0; t--){
        single::solve_single_to_single(t,sol,par); 
        single::solve_couple_to_single(t,sol,par); 
        couple::solve_couple(t,sol,par);
        couple::solve_single_to_couple(t,sol,par);
        single::expected_value_start_single(t,sol,par);
    }
}


EXPORT void simulate(sim_struct *sim, sol_struct *sol, par_struct *par){
    
    sim::model(sim,sol,par);

}


EXPORT void compute_margEV(sol_struct* sol, par_struct* par){
    for (int t = 0; t < par->T; t++){
        single::calc_marginal_value_single(t, woman, sol, par);
        single::calc_marginal_value_single(t, man, sol, par);

        for (int iP=0; iP<par->num_power; iP++){
            for (int iL=0; iL<par->num_love; iL++){
                int idx = index::couple(t,iP,iL,0,par);
                double* EVw = &sol->EVw_start_as_couple[idx];
                double* EVm = &sol->EVm_start_as_couple[idx];
                double* EmargV = &sol->EmargV_start_as_couple[idx];
                couple::calc_marginal_value_couple(t, iP, iL, EVw, EVm, EmargV, sol, par);
            }
        }
    }
}


EXPORT double calc_init_mu(int t, double love, double Aw, double Am, sol_struct* sol, par_struct* par){
    logs::write("barg_log.txt", 0, "calc_init_mu\n");
    double power =  single::calc_initial_bargaining_weight(t, love, Aw, Am, sol, par);
    logs::write("barg_log.txt", 1, "poewr: %f\n", power);
    return power;
}