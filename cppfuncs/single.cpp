#ifndef MAIN
#define SINGLE
#include "myheader.cpp"
#endif

namespace single {
    typedef struct {
        
        double M;
        double *V_next;
        int gender;
        par_struct *par;

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

    double objfunc_single(unsigned n, const double *x, double *grad, void *solver_data_in){
        double love = 0.0;

        // unpack
        solver_single_struct *solver_data = (solver_single_struct *) solver_data_in;
        
        double C_tot = x[0];
        int gender = solver_data->gender;
        double M = solver_data->M;
        par_struct *par = solver_data->par;

        // flow-utility
        double C_priv = cons_priv_single(C_tot,gender,par);
        double C_pub = C_tot - C_priv;
        
        double Util = utils::util(C_priv,C_pub,gender,par,love);
        
        // continuation value
        double *grid_A = par->grid_Aw;
        if (gender==man){
            grid_A = par->grid_Am;
        }
        double A = M - C_tot;

        double Vnext = tools::interp_1d(grid_A,par->num_A,solver_data->V_next,A);
        
        // return discounted sum
        return -(Util + par->beta*Vnext);
    }

    void solve_single(int t,sol_struct *sol,par_struct *par){
        double love = 0.0; // no love for singles

        // terminal period
        if (t == (par->T-1)){
            for (int iA=0; iA<par->num_A;iA++){
                int idx = index::index2(t,iA,par->T,par->num_A);

                double Aw = par->grid_Aw[iA];
                double Am = par->grid_Am[iA];

                double Cw = par->R*Aw + par->inc_w;
                double Cm = par->R*Am + par->inc_m;
                
                sol->Cw_priv_single[idx] = cons_priv_single(Cw,woman,par);
                sol->Cw_pub_single[idx] = Cw - sol->Cw_priv_single[idx];
                sol->Vw_single[idx] = utils::util(sol->Cw_priv_single[idx],sol->Cw_pub_single[idx],woman,par,love);
                
                sol->Cm_priv_single[idx] = cons_priv_single(Cm,man,par);
                sol->Cm_pub_single[idx] = Cm - sol->Cm_priv_single[idx];
                sol->Vm_single[idx] = utils::util(sol->Cm_priv_single[idx],sol->Cm_pub_single[idx],man,par,love);
            }
        } else {
            
            #pragma omp parallel num_threads(par->threads)
            {

                // 1. allocate objects for solver
                solver_single_struct* solver_data = new solver_single_struct;
                
                int dim = 1;
                double lb[1],ub[1],x[1];
                
                auto opt = nlopt_create(NLOPT_LN_BOBYQA, dim); // NLOPT_LD_MMA NLOPT_LD_LBFGS NLOPT_GN_ORIG_DIRECT
                double minf=0.0;

                // 2. loop over assets
                #pragma omp for
                for (int iA=0; iA<par->num_A;iA++){
                    int idx = index::index2(t,iA,par->T,par->num_A);

                    double Aw = par->grid_Aw[iA];
                    double Am = par->grid_Am[iA];
                    
                    // resources
                    double Mw = par->R*Aw + par->inc_w;
                    double Mm = par->R*Am + par->inc_m;
                    
                    // search over optimal total consumption, C
                    // WOMEN
                    // settings
                    solver_data->M = Mw;
                    solver_data->V_next = &sol->Vw_single[index::index2(t+1,0,par->T,par->num_A)];
                    solver_data->gender = woman;
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

                    // store results
                    double Cw = x[0];
                    sol->Cw_priv_single[idx] = cons_priv_single(Cw,woman,par);
                    sol->Cw_pub_single[idx] = Cw - sol->Cw_priv_single[idx];
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
                    sol->Cm_priv_single[idx] = cons_priv_single(Cm,man,par);
                    sol->Cm_pub_single[idx] = Cm - sol->Cm_priv_single[idx];
                    sol->Vm_single[idx] = -minf;            
                    
                } // iA

                // 4. destoy optimizer
                nlopt_destroy(opt);

            } // pragma
        }   
        
    }


}