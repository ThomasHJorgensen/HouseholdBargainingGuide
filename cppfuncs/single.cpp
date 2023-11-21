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
        double constant = alpha1*pow(share,phi) + alpha2*pow(1.0-share,phi);
        return phi * pow(C_tot,(1.0-rho)*phi -1.0 ) * pow(constant,1.0 - rho);

    }
    double resources(double A,int gender,par_struct* par) {
        double income = par->inc_w;
        if (gender == man) {
            income = par->inc_m;
        }
        return par->R*A + income;
    }

    double value_of_choice(double C_tot,double M, int gender, double* V_next, par_struct* par){

        // flow-utility
        double Util = util_C(C_tot,gender,par);
        // logs::write("egm_singles.txt",1,"Util = %f, gender = %f\n",Util,gender);
        
        // continuation value
        double *grid_A = par->grid_Aw; 
        if (gender==man){
            grid_A = par->grid_Am;
        }
        double A = M - C_tot;

        double Vnext = tools::interp_1d(grid_A,par->num_A,V_next,A);
        
        // return discounted sum
        return Util + par->beta*Vnext;
    }

    double objfunc_single(unsigned n, const double *x, double *grad, void *solver_data_in){
        double love = 0.0;

        // unpack
        solver_single_struct *solver_data = (solver_single_struct *) solver_data_in;  
        
        double C_tot = x[0];
        int gender = solver_data->gender;
        double M = solver_data->M;
        par_struct *par = solver_data->par;

        return - value_of_choice(C_tot,M,gender,solver_data->V_next,par);

    }

    void EGM_step_single(int t, sol_struct *sol, par_struct *par){
        for (int iA=0; iA<par->num_A_pd; iA++){
            logs::write("egm_singles.txt",1,"t = %d, iA = %d\n",t,iA);

            // Unpack
            int idx      = index::index2(t  ,iA,par->T,par->num_A);
            int idx_next = index::index2(t+1,iA,par->T,par->num_A);
            double love = 0.0;
            const double A = par->grid_A_pd[iA];


            // // Marginal value of next period wealth
            // double w = par->beta * sol->

            // Next period consumption - Note: can change index to just t+1
            double C_next_w = tools::interp_1d(par->grid_Aw, par->num_A,&sol->Cw_tot_single[index::index2(t+1,0,par->T,par->num_A)], A);
            double C_next_m = tools::interp_1d(par->grid_Am, par->num_A,&sol->Cm_tot_single[index::index2(t+1,0,par->T,par->num_A)], A);

            // logs::write("egm_singles.txt",1,"C_next_w = %f, C_next_m = %f\n",C_next_w,C_next_m);

            // Marginal utility of next period consumption
            // I need expectation around here
            double EmargUw = par->beta * par->R * marg_util_C(C_next_w, woman, par);
            double EmargUm = par->beta * par->R * marg_util_C(C_next_m, man,   par);
            // logs::write("egm_singles.txt",1,"EmargUw = %f, EmargUm = %f\n",EmargUw,EmargUm);

            // // alternative
            // double EmargUw = sol->EmargUw_single_pd[index::index2(t+1,iA,par->T,par->num_A)];
            // double EmargUm = sol->EmargUm_single_pd[index::index2(t+1,iA,par->T,par->num_A)];

            // Consumption today
            double C_now_w = tools::interp_1d(par->grid_marg_u_single_w_for_inv, par->num_Ctot, par->grid_inv_marg_u, EmargUw);
            double C_now_m = tools::interp_1d(par->grid_marg_u_single_m_for_inv, par->num_Ctot, par->grid_inv_marg_u, EmargUm);
            // logs::write("egm_singles.txt",1,"C_now_w = %f, C_now_m = %f\n",C_now_w,C_now_m);

            // Asset grid today
            double M_now_w = A + C_now_w;
            double M_now_m = A + C_now_m;

            // Store
            sol->C_totw_single_pd[iA] = C_now_w;
            sol->C_totm_single_pd[iA] = C_now_m;

            sol->Mw_single_pd[iA] = M_now_w;
            sol->Mm_single_pd[iA] = M_now_m;
            // logs::write("egm_singles.txt",1,"M_now_w = %f, M_now_m = %f\n",A_now_w,A_now_m);


        }

        for (int iA=0; iA<par->num_A; iA++){
            int idx = index::index2(t,iA,par->T,par->num_A);
            int idx_next = index::index2(t+1,iA,par->T,par->num_A);

            // resources
            double Mw = resources(par->grid_Aw[iA], woman, par);
            double Mm = resources(par->grid_Am[iA],   man, par);
            logs::write("egm_singles.txt",1,"Mw = %f, Mm = %f\n",Mw,Mm);

            // interp back to exogenous grid
            sol->Cw_tot_single[idx] = tools::interp_1d(sol->Mw_single_pd, par->num_A_pd, sol->C_totw_single_pd, Mw);
            sol->Cm_tot_single[idx] = tools::interp_1d(sol->Mm_single_pd, par->num_A_pd, sol->C_totm_single_pd, Mm);
            // logs::write("egm_singles.txt",1,"Cw_tot_single = %f, Cm_tot_single = %f\n",sol->C_totw_single_pd[idx],sol->C_totm_single_pd[idx]);
            // logs::write("egm_singles.txt",1,"C_totw_single_pd = %f, C_totm_single_pd = %f\n",sol->C_totw_single_pd[idx],sol->C_totm_single_pd[idx]);
            // logs::write("egm_singles.txt",1,"Cw_tot_single = %f, Cm_tot_single = %f\n",sol->Cw_tot_single[idx],sol->Cm_tot_single[idx]);

            // Split consumption
            intraperiod_allocation(&sol->Cw_priv_single[idx],&sol->Cw_pub_single[idx],sol->Cw_tot_single[idx],woman,par);
            intraperiod_allocation(&sol->Cm_priv_single[idx],&sol->Cm_pub_single[idx],sol->Cm_tot_single[idx],  man,par);

            // values
            sol->Vw_single[idx] = value_of_choice(sol->Cw_tot_single[idx], Mw, woman, &sol->Vw_single[index::index2(t+1,0,par->T,par->num_A)], par);
            sol->Vm_single[idx] = value_of_choice(sol->Cm_tot_single[idx], Mm,   man, &sol->Vm_single[index::index2(t+1,0,par->T,par->num_A)], par);
        }
    }


    void solve_single(int t,sol_struct *sol,par_struct *par){
        double love = 0.0; // no love for singles 

        // terminal period
        if (t == (par->T-1)){
            logs::write("egm_singles.txt",0,"");

            for (int iA=0; iA<par->num_A;iA++){
                int idx = index::index2(t,iA,par->T,par->num_A);

                double Aw = par->grid_Aw[iA];
                double Am = par->grid_Am[iA];

                double Cw = resources(Aw,woman,par); 
                double Cm = resources(Am,man,par); 

                // Save consumption
                sol->Cw_tot_single[idx] = Cw;
                sol->Cm_tot_single[idx] = Cm;
                
                intraperiod_allocation(&sol->Cw_priv_single[idx],&sol->Cw_pub_single[idx],Cw,woman,par);
                sol->Vw_single[idx] = utils::util(sol->Cw_priv_single[idx],sol->Cw_pub_single[idx],woman,par,love);
                
                intraperiod_allocation(&sol->Cm_priv_single[idx],&sol->Cm_pub_single[idx],Cm,man,par);
                sol->Vm_single[idx] = utils::util(sol->Cm_priv_single[idx],sol->Cm_pub_single[idx],man,par,love);
            }
        } else {
            
            #pragma omp parallel num_threads(par->threads)
            {
                if (par->do_egm){
                    EGM_step_single(t, sol, par);
                } else {
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
                } // else
            } // pragma
        }   
        
    }


}