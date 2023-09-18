#ifndef MAIN
#define SINGLE
#include "myheader.cpp"
#endif

namespace single {
    typedef struct { //AMO: define structure to hold data relevant for optimizer
        
        double M;               //AMO: resources
        double *V_next;         //AMO: pointer to continuation value
        int gender;             //AMO: gender indicator
        par_struct *par;        //AMO: pointer to par
                                //AMO Q: not time? sol?

    } solver_single_struct;


    double cons_priv_single(double C_tot,int gender,par_struct *par){
        // closed form solution for intra-period problem of single.
        double rho = par->rho_w;
        double phi = par->phi_w;
        double alpha1 = par->alpha1_w;
        double alpha2 = par->alpha2_w;
        if (gender == man) {
            rho = par->rho_m;
            phi = par->phi_m;
            alpha1 = par->alpha1_m;
            alpha2 = par->alpha2_m;
        }  
        
        return C_tot/(1.0 + pow(alpha2/alpha1,1.0/(1.0-phi) ));
    }

    void intraperiod_allocation(double* C_priv, double* C_pub , double C_tot, int gender,par_struct *par){
        C_priv[0] = cons_priv_single(C_tot,gender,par);
        C_pub[0] = C_tot - C_priv[0];
    }

    double util_C(double C_tot, int gender, par_struct* par){
        double love = 0.0;
        
        // flow-utility
        double C_priv = cons_priv_single(C_tot,gender,par);
        double C_pub = C_tot - C_priv;
        
        return utils::util(C_priv,C_pub,gender,par,love);
    }

    double marg_util_C(double C_tot, int gender, par_struct* par){
        //AMO: closed for marginal utility of C_tot
        double rho = par->rho_w;
        double phi = par->phi_w;
        double alpha1 = par->alpha1_w;
        double alpha2 = par->alpha2_w;
        if (gender == man) {
            rho = par->rho_m;
            phi = par->phi_m;
            alpha1 = par->alpha1_m;
            alpha2 = par->alpha2_m;
        }  
        
        double share = 1.0/(1.0 + pow(alpha2/alpha1,1.0/(1.0-phi) ));
        double constant = alpha1*pow(share,phi) + alpha2*pow(1.0-constant,phi);
        return phi * pow(C_tot,(1.0-rho)*phi -1.0 ) * pow(constant,1.0 - rho);

    }
    double resources(double A,int gender,par_struct* par) {
        //AMO: resources when single
        double income = par->inc_w;
        if (gender == man) {
            income = par->inc_m;
        }
        return par->R*A + income;
    }

    double value_of_choice(double C_tot,double M, int gender, double* V_next, par_struct* par){

        // flow-utility
        double Util = util_C(C_tot,gender,par);
        
        // continuation value
        double *grid_A = par->grid_Aw; //AMO: pointer to relevant assets grid (w/m)
        if (gender==man){
            grid_A = par->grid_Am;
        }
        double A = M - C_tot;

        double Vnext = tools::interp_1d(grid_A,par->num_A,V_next,A);
        
        // return discounted sum
        return Util + par->beta*Vnext;
    }

    double objfunc_single(unsigned n, const double *x, double *grad, void *solver_data_in){
                //AMO: objective function takes (const double &x, double &grad, void* solver_data) as input according to NLopt docs
                // - and appearantly also n (dimension)
                //AMO: unsigned means only non-negative values (integers)
                //AMO: const keyword indicates that x should not be modified by this function
                //AMO Q: grad is not used - just for optimizer syntax?
        double love = 0.0;

        // unpack
        solver_single_struct *solver_data = (solver_single_struct *) solver_data_in;  //AMO: casts solver_data_in to solver_data (which is type solver_single_struct)
        
        double C_tot = x[0];
        int gender = solver_data->gender;
        double M = solver_data->M;
        par_struct *par = solver_data->par;

        return - value_of_choice(C_tot,M,gender,solver_data->V_next,par);

    }

    void solve_single(int t,sol_struct *sol,par_struct *par){
        double love = 0.0; // no love for singles  //AMO: ouch

        // terminal period
        if (t == (par->T-1)){
            for (int iA=0; iA<par->num_A;iA++){
                int idx = index::index2(t,iA,par->T,par->num_A);

                double Aw = par->grid_Aw[iA];
                double Am = par->grid_Am[iA];

                double Cw = resources(Aw,woman,par); //AMO: consume everything in final period
                double Cm = resources(Am,man,par); 
                
                intraperiod_allocation(&sol->Cw_priv_single[idx],&sol->Cw_pub_single[idx],Cw,woman,par);
                sol->Vw_single[idx] = utils::util(sol->Cw_priv_single[idx],sol->Cw_pub_single[idx],woman,par,love);
                
                intraperiod_allocation(&sol->Cm_priv_single[idx],&sol->Cm_pub_single[idx],Cm,man,par);
                sol->Vm_single[idx] = utils::util(sol->Cm_priv_single[idx],sol->Cm_pub_single[idx],man,par,love);
            }
        } else {
            
            #pragma omp parallel num_threads(par->threads) //AMO: settings for parallelization
            {

                // 1. allocate objects for solver
                solver_single_struct* solver_data = new solver_single_struct;
                
                int dim = 1;
                double lb[1],ub[1],x[1];
                
                auto opt = nlopt_create(NLOPT_LN_BOBYQA, dim); // NLOPT_LD_MMA NLOPT_LD_LBFGS NLOPT_GN_ORIG_DIRECT
                // AMO: opt objects takes (optimization algorithm, dimensions) as input
                double minf=0.0; //AMO: initialize: will eventually contain the value of the objective fct at minimum

                // 2. loop over assets
                #pragma omp for
                for (int iA=0; iA<par->num_A;iA++){
                    int idx = index::index2(t,iA,par->T,par->num_A);
                    
                    // resources
                    double Aw = par->grid_Aw[iA];
                    double Am = par->grid_Am[iA];
                    
                    double Mw = resources(Aw,woman,par); 
                    double Mm = resources(Am,man,par); 
                    
                    // search over optimal total consumption, C
                    // WOMEN
                    // settings
                    solver_data->M = Mw;
                    solver_data->V_next = &sol->Vw_single[index::index2(t+1,0,par->T,par->num_A)];
                    solver_data->gender = woman;
                    solver_data->par = par; //AMO: update solver_data
                    nlopt_set_min_objective(opt, objfunc_single, solver_data); 
                    //AMO: set objective function: find min of objfunc_single with solver_data using opt optimizer
                        
                    // bounds
                    lb[0] = 1.0e-8;
                    ub[0] = solver_data->M;
                    nlopt_set_lower_bounds(opt, lb);
                    nlopt_set_upper_bounds(opt, ub);

                    // optimize
                    x[0] = solver_data->M/2.0; //AMO: set starting value 
                    nlopt_optimize(opt, x, &minf); //AMO: do optimization. Takes (function, &starting val, &fct val)

                    // store results
                    double Cw = x[0]; //AMO: optimizer result
                    //AMO Q: Is there a way to verify succesfull optimization?
                    intraperiod_allocation(&sol->Cw_priv_single[idx],&sol->Cw_pub_single[idx],Cw,woman,par);
                    sol->Vw_single[idx] = -minf;

                    // MEN
                    // settings
                    solver_data->M = Mm;
                    solver_data->V_next = &sol->Vm_single[index::index2(t+1,0,par->T,par->num_A)];
                    solver_data->gender = man;
                    solver_data->par = par;
                    nlopt_set_min_objective(opt, objfunc_single, solver_data);
                        
                    // bounds
                    lb[0] = 1.0e-8;
                    ub[0] = solver_data->M;
                    nlopt_set_lower_bounds(opt, lb);
                    nlopt_set_upper_bounds(opt, ub);

                    // optimize
                    x[0] = solver_data->M/2.0;
                    nlopt_optimize(opt, x, &minf);

                    double Cm = x[0];
                    intraperiod_allocation(&sol->Cm_priv_single[idx],&sol->Cm_pub_single[idx],Cm,man,par);
                    sol->Vm_single[idx] = -minf;            
                    
                } // iA

                // 4. destroy optimizer
                nlopt_destroy(opt);

            } // pragma
        }   
        
    }


}