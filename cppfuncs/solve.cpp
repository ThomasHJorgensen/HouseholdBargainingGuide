#define MAIN
#include "myheader.h"


typedef struct { //AMO: namespace type structure
    double power;
    double C_tot;

    par_struct *par;

} solver_precompute_struct;

double objfunc_precompute(unsigned n, const double *x, double *grad, void *solver_data_in){
    // unpack
    solver_precompute_struct *solver_data = (solver_precompute_struct *) solver_data_in; //AMO: casts (pointer to) solver_data_in to (pointer to) solver_data (of type solver_precompute_struct)
    
    double C_tot = solver_data->C_tot;
    double power = solver_data->power;
    par_struct *par = solver_data->par;

    double Cw_priv = x[0];  // AMO: bc x is pointer 
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
    solver_precompute_struct* solver_data = new solver_precompute_struct;  // AMO: allocates memory for new instance of solver_precompute_struct
            
    int dim = 2;
    double lb[2],ub[2],x[2];   // AMO: arrays of size 2, settings for optimizer
    
    auto opt = nlopt_create(NLOPT_LN_BOBYQA, dim); // NLOPT_LD_MMA NLOPT_LD_LBFGS NLOPT_GN_ORIG_DIRECT
    double minf=0.0;

    // search over optimal total consumption, C
    // settings
    solver_data->C_tot = C_tot;         //AMO: fill in solver_data
    solver_data->power = power;
    solver_data->par = par;
    nlopt_set_min_objective(opt, objfunc_precompute, solver_data);   //AMO: pass objective and data to optimizer object (see EconModel guide)
    nlopt_set_maxeval(opt, 2000);

    // bounds
    lb[0] = 0.0;                //AMO: optimizer settings
    lb[1] = 0.0;
    ub[0] = solver_data->C_tot;
    ub[1] = solver_data->C_tot;
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);

    // optimize
    x[0] = solver_data->C_tot/3.0;
    x[1] = solver_data->C_tot/3.0;
    nlopt_optimize(opt, x, &minf);          //AMO: run optimizer to find optimal intra-temporal allocation
    nlopt_destroy(opt);                 //AMO Q: why adress of minf but just x? 
    
    // unpack
    Cw_priv[0] = x[0];
    Cm_priv[0] = x[1];
    C_pub[0] = C_tot - Cw_priv[0] - Cm_priv[0];

}

void precompute(sol_struct* sol, par_struct* par){
    #pragma omp parallel num_threads(par->threads)      //AMO: parallelization settings
    {
        #pragma omp for         //AMO: do parallelized loop 
        for (int iP=0; iP<par->num_power; iP++){  //AMO: for iP in range(num_power)
            for (int i=0; i < par->num_Ctot; i++){
                double C_tot = par->grid_Ctot[i];
                int idx = index::index2(iP,i,par->num_power,par->num_Ctot);         //AMO: index for iP, i in power X Ctot grid
                solve_intraperiod_couple(&sol->pre_Ctot_Cw_priv[idx], &sol->pre_Ctot_Cm_priv[idx], &sol->pre_Ctot_C_pub[idx] , C_tot,par->grid_power[iP],par);
                                    // AMO: &sol->xx means return the adress of sol->xx (AMO Q: should & always be used when the fct input is defined as pointer?)
                
                // calculate marginal utility and inverse marginal utility for EGM
                if(par->do_egm){
                    int iL = 0; // does not matter for the marginal utility //AMO: love index 
                    int idx_lag = index::index2(iP,i-1,par->num_power,par->num_Ctot);       // AMO: index for consumption at point below (for finite dif approx of marg U)

                    // utility at current allocation
                    par->grid_util[idx] = utils::util_couple(sol->pre_Ctot_Cw_priv[idx],sol->pre_Ctot_Cm_priv[idx],sol->pre_Ctot_C_pub[idx],iP,iL,par);
                    
                    if (i>0) {

                        // marginal utility. Use finite difference but closed form could be used.
                        par->grid_marg_u[idx_lag] = (par->grid_util[idx] - par->grid_util[idx_lag])/(C_tot - par->grid_Ctot[i-1]);
                     
                        // inverse marginal utility: flip the grid of marginal util (such that ascending) and store as new "x-axis" grid
                        //AMO: this is necessary for the binary search in the interpolaiton algo
                        int idx_flip_lag = index::index2(iP,par->num_Ctot-1 - (i-1),par->num_power,par->num_Ctot);
                        par->grid_marg_u_for_inv[idx_flip_lag] = par->grid_marg_u[idx_lag];

                        if (i==(par->num_Ctot-1)){
                            par->grid_marg_u[idx] = par->grid_marg_u[idx_lag]; // impose constant slope at last grid-point

                            int idx_flip = index::index2(iP,par->num_Ctot-1 - i,par->num_power,par->num_Ctot);
                            par->grid_marg_u_for_inv[idx_flip] = par->grid_marg_u[idx_lag];

                        }

                    }

                } // EGM 

            }
        }
    }
}

/////////////
// 5. MAIN //
/////////////

EXPORT void solve(sol_struct *sol, par_struct *par){
    
    // pre-compute intra-temporal optimal allocation
    precompute(sol,par);

    // loop backwards
    for (int t = par->T-1; t >= 0; t--){

        single::solve_single(t,sol,par); //AMO: run solve_single fct from namespace single, uses backwards induction
        couple::solve_couple(t,sol,par);

    }
}


EXPORT void simulate(sim_struct *sim, sol_struct *sol, par_struct *par){
    
    sim::model(sim,sol,par);

}
