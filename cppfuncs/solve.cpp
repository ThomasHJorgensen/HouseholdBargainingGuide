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
    }
}


EXPORT void simulate(sim_struct *sim, sol_struct *sol, par_struct *par){
    
    sim::model(sim,sol,par);

}



EXPORT double draw_partner_assets(double A, int gender, int i, int t, sim_struct *sim, par_struct *par){
    return sim::draw_partner_assets(A,gender,i,t,sim,par);
}