#define MAIN
#include "myheader.h"


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

        // if(t>=(par->T-1)){
        //     single::solve_single_to_single(t,sol,par);
        //     single::solve_couple_to_single(t,sol,par); //AMO: maybe put this in couple namespace?
        //     couple::solve_remain_start_trans_couple(t,sol,par); 
        //     single::expected_value_start_single(t,sol,par); // expectation over potential partners, using value of trans to couple and value of remain single
        // }
        // if(t>=(par->T-2)){
        //     single::solve_remain_trans_single(t,sol,par); // TODO: introduce remain_single solution and update code.
        //     couple::solve_remain_start_trans_couple(t,sol,par); 
        //     // single::expected_value_start_single(t,sol,par); // expectation over potential partners, using value of trans to couple and value of remain single
        // }

    }
}


EXPORT void simulate(sim_struct *sim, sol_struct *sol, par_struct *par){
    
    sim::model(sim,sol,par);

}



