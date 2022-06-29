#define MAIN
#include "myheader.h"


typedef struct {
    double power;
    double C_tot;

    par_struct *par;

} solver_precompute_struct;

double objfunc_precompute(unsigned n, const double *x, double *grad, void *solver_data_in){
    // unpack
    solver_precompute_struct *solver_data = (solver_precompute_struct *) solver_data_in;
    
    double C_tot = solver_data->C_tot;
    double power = solver_data->power;
    par_struct *par = solver_data->par;

    double Cw_priv = x[0];
    double Cm_priv = x[1];
    double C_pub = C_tot - Cw_priv - Cm_priv;

    // weighted utility of choice
    double love = 0.0; // does not matter for optimal allocation
    double val = power*utils::util(Cw_priv,C_pub,woman,par,love) + (1.0-power)*utils::util(Cm_priv,C_pub,man,par,love);

    // return negative of value
    return - val;
}

void solve_intraperiod_couple(double* Cw_priv,double* Cm_priv,double* C_pub , double C_tot,double power,par_struct *par){
    
    // setup numerical solver
    solver_precompute_struct* solver_data = new solver_precompute_struct;
            
    int dim = 2;
    double lb[2],ub[2],x[2];
    
    auto opt = nlopt_create(NLOPT_LN_BOBYQA, dim); // NLOPT_LD_MMA NLOPT_LD_LBFGS NLOPT_GN_ORIG_DIRECT
    double minf=0.0;

    // search over optimal total consumption, C
    // settings
    solver_data->C_tot = C_tot;
    solver_data->power = power;
    solver_data->par = par;
    nlopt_set_min_objective(opt, objfunc_precompute, solver_data);
    nlopt_set_maxeval(opt, 2000);

    // bounds
    lb[0] = 0.0;
    lb[1] = 0.0;
    ub[0] = solver_data->C_tot;
    ub[1] = solver_data->C_tot;
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);

    // optimize
    x[0] = solver_data->C_tot/3.0;
    x[1] = solver_data->C_tot/3.0;
    nlopt_optimize(opt, x, &minf);
    nlopt_destroy(opt);
    
    // unpack
    Cw_priv[0] = x[0];
    Cm_priv[0] = x[1];
    C_pub[0] = C_tot - Cw_priv[0] - Cm_priv[0];

}

/////////////
// 5. MAIN //
/////////////

EXPORT void solve(sol_struct *sol, par_struct *par){
    
    // pre-compute intra-temporal optimal allocation
    #pragma omp parallel num_threads(par->threads)
    {
        #pragma omp for
        for (int iP=0; iP<par->num_power; iP++){
            for (int i=0; i < par->num_Ctot; i++){
                int idx = index::index2(iP,i,par->num_power,par->num_Ctot);
                solve_intraperiod_couple(&sol->pre_Ctot_Cw_priv[idx], &sol->pre_Ctot_Cm_priv[idx], &sol->pre_Ctot_C_pub[idx] , par->grid_Ctot[i],par->grid_power[iP],par);
            }
        }
    }

    // loop backwards
    for (int t = par->T-1; t >= 0; t--){

        single::solve_single(t,sol,par);
        couple::solve_couple(t,sol,par);

    }
}


EXPORT void simulate(sim_struct *sim, sol_struct *sol, par_struct *par){
    
    sim::model(sim,sol,par);

}
