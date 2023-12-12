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

    void intraperiod_allocation(double* C_priv, double* C_pub , double C_tot, int gender,par_struct *par){
        C_priv[0] = utils::cons_priv_single(C_tot,gender,par);
        C_pub[0] = C_tot - C_priv[0];
    }

    // }
    double resources(double A,int gender,par_struct* par) {
        double income = par->inc_w;
        if (gender == man) {
            income = par->inc_m;
        }
        return par->R*A + income;
    }

    double value_of_choice(double C_tot,double M, int gender, double* V_next, par_struct* par){

        // flow-utility
        double Util = utils::util_C_single(C_tot,gender,par);
        
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

    void EGM_single(int t, int gender, sol_struct* sol, par_struct* par){
        // 1. Setup
        /// a. unpack
        int const &T {par->T};
        int const &num_marg_u {par->num_marg_u};
        int const &num_A {par->num_A};
        int const &num_A_pd {par->num_A_pd};
        bool const &analytic_inv_marg_u_single {par->analytic_inv_marg_u_single};
        double const &beta {par->beta};
        double const &R {par->R};
        double* const &grid_inv_marg_u {par->grid_inv_marg_u};

        /// b. gender specific variables
        //// o. woman
        double* grid_A {par->grid_Aw};
        double* grid_A_pd {par->grid_Aw_pd};
        double* grid_marg_u_single_for_inv {par->grid_marg_u_single_w_for_inv};
        double* V {sol->Vw_single};
        double* marg_V {sol->marg_Vw_single};
        double* C_tot {sol->Cw_tot_single};
        double* C_priv {sol->Cw_priv_single};
        double* C_pub {sol->Cw_pub_single};
        //// oo. man
        if (gender == man){
            grid_A = par->grid_Am;
            grid_A_pd = par->grid_Am_pd;
            grid_marg_u_single_for_inv = par->grid_marg_u_single_m_for_inv;
            V = sol->Vm_single;
            marg_V = sol->marg_Vm_single;
            C_tot = sol->Cm_tot_single;
            C_priv = sol->Cm_priv_single;
            C_pub = sol->Cm_pub_single;
        }

        /// c. Allocate memory
        double* EmargU_pd {new double[num_A_pd]};
        double* C_tot_pd {new double[num_A_pd]};
        double* M_pd {new double[num_A_pd]};

        // 2. EGM step
        /// setup
        int idx_next = index::index2(t+1,0, T, num_A);
        int min_point_A {0};

        for (int iA_pd=0; iA_pd<num_A_pd; iA_pd++){

            /// a. get next period assets
            double &A_next = grid_A_pd[iA_pd];

            /// b. calculate expected marginal utility
            min_point_A = tools::binary_search(min_point_A, num_A, grid_A, A_next);
            EmargU_pd[iA_pd] = tools::interp_1d_index(grid_A, num_A, &marg_V[idx_next],A_next, min_point_A);

            /// c. invert marginal utility by interpolation from pre-computed grid
            if (analytic_inv_marg_u_single == 1){
                C_tot_pd[iA_pd] = utils::inv_marg_util_C(EmargU_pd[iA_pd], gender, par);
            } else {
                C_tot_pd[iA_pd] = tools::interp_1d(grid_marg_u_single_for_inv, num_marg_u, grid_inv_marg_u, EmargU_pd[iA_pd]);
            }
            /// d. endogenous grid over resources
            M_pd[iA_pd] = C_tot_pd[iA_pd] + A_next;
        }

        // 3. interpolate to common grid
        ///setup
        int idx = index::index2(t,0, T, num_A);
        min_point_A = 0;

        for (int iA=0; iA<num_A; iA++){

            /// a. calculate resources
            double M_now = resources(grid_A[iA], gender, par);

            /// b. find total consumption
            min_point_A = tools::binary_search(min_point_A,num_A_pd, M_pd, M_now);
            C_tot[idx] = tools::interp_1d_index(M_pd, num_A_pd, C_tot_pd, M_now, min_point_A);

            /// c. handle credit constraint 
            //// if credit constrained
            if (M_now <= M_pd[0]) {
                ///// o. consume all resources
                C_tot[idx] = M_now; 

                ///// oo. calculate marginal value of constrained consumption
                if (par->analytic_marg_u_single){
                    marg_V[idx] = beta * R * utils::marg_util_C(C_tot[idx], gender, par);
                }
                else{
                    marg_V[idx] = beta * R * tools::interp_1d(par->grid_C_for_marg_u, par->num_marg_u, par->grid_marg_u, C_tot[idx]);
                }
            }
            //// if not credit constrained
            else{
                // o. calculate marginal value of unconstrained consumption
                marg_V[idx] = beta * R * tools::interp_1d_index(M_pd, num_A_pd, EmargU_pd, M_now, min_point_A);
            }

            /// d. calculate private and public consumption
            intraperiod_allocation(&C_priv[idx], &C_pub[idx], C_tot[idx], gender, par);
            
            /// e. calculate values (not used in EGM step)
            V[idx] = value_of_choice(C_tot[idx], M_now, gender, &V[idx_next], par);

            /// f. update index
            idx++;
        }

        // 4. clean up
        delete[] EmargU_pd;
        delete[] C_tot_pd;
        delete[] M_pd;
    }



    void solve_single(int t,sol_struct *sol,par_struct *par){
        double love = 0.0; // no love for singles 

        // terminal period
        if (t == (par->T-1)){
            for (int iA=0; iA<par->num_A;iA++){
                int idx = index::index2(t,iA,par->T,par->num_A);

                double Aw = par->grid_Aw[iA];
                double Am = par->grid_Am[iA];

                sol->Cw_tot_single[idx] = resources(Aw,woman,par); 
                sol->Cm_tot_single[idx] = resources(Am,man,par); 
                
                intraperiod_allocation(&sol->Cw_priv_single[idx],&sol->Cw_pub_single[idx],sol->Cw_tot_single[idx],woman,par);
                sol->Vw_single[idx] = utils::util(sol->Cw_priv_single[idx],sol->Cw_pub_single[idx],woman,par,love);
                
                intraperiod_allocation(&sol->Cm_priv_single[idx],&sol->Cm_pub_single[idx],sol->Cm_tot_single[idx],man,par);
                sol->Vm_single[idx] = utils::util(sol->Cm_priv_single[idx],sol->Cm_pub_single[idx],man,par,love);

                if (par->do_egm) {
                    sol->marg_Vw_single[idx] = par->beta*par->R*utils::marg_util_C(sol->Cw_tot_single[idx], woman, par);
                    sol->marg_Vm_single[idx] = par->beta*par->R*utils::marg_util_C(sol->Cm_tot_single[idx], man, par);
                }
            }
        } else {
            if (par->do_egm) {
                EGM_single(t, woman, sol, par);
                EGM_single(t, man, sol, par);
            }
            else {
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
                        sol->Cw_tot_single[idx] = x[0];
                        intraperiod_allocation(&sol->Cw_priv_single[idx],&sol->Cw_pub_single[idx],sol->Cw_tot_single[idx],woman,par);
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

                        sol->Cm_tot_single[idx] = x[0];
                        intraperiod_allocation(&sol->Cm_priv_single[idx],&sol->Cm_pub_single[idx],sol->Cm_tot_single[idx],man,par);
                        sol->Vm_single[idx] = -minf;       
                        
                    } // iA

                    // 4. destroy optimizer
                    nlopt_destroy(opt);

                } // pragma
            }
        }   
        
    }


}