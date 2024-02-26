// functions for bargaining process.
#ifndef MAIN
#define BARGAINING
#include "myheader.cpp"
#endif

namespace bargaining {

    int find_left(double *S, int Nx){

        int left_point {};

        if (S[0] <= S[Nx-1]){
            left_point = tools::binary_search(0, Nx, S, 0.0);
        }
        else{
            left_point = tools::binary_search_over_descending_function(0, Nx, S, 0.0);
        }

        return left_point;
}

    // divorce
    void divorce(int iP, int *power_idx, double *power, index::index_couple_struct *idx_couple, double **list_start_as_couple, double *list_trans_to_single, int num, par_struct *par){
        int idx = idx_couple->idx(iP);
        power_idx[idx] = -1;
        power[idx] = -1.0;

        for (int i = 0; i < num; i++){
            list_start_as_couple[i][idx] = list_trans_to_single[i];
        }
    }
    
    // remain
    void remain(int iP, int *power_idx, double *power, index::index_couple_struct *idx_couple, double **list_start_as_couple, double **list_remain_couple, int num, par_struct *par){
        int idx = idx_couple->idx(iP);
        power_idx[idx] = iP;
        power[idx] = par->grid_power[iP];

        for (int i = 0; i < num; i++){
            list_start_as_couple[i][idx] = list_remain_couple[i][idx];
        }
    }


    // update to indifference point
    void update_to_indifference(int iP, int left_point, int low_point, double power_at_zero, int *power_idx, double *power, index::index_couple_struct *idx_couple, double **list_start_as_couple, double **list_remain_couple, int num, par_struct *par, int sol_idx = -1){
        int idx = idx_couple->idx(iP);
        power_idx[idx] = low_point;
        power[idx] = power_at_zero;

        int delta = idx_couple->idx(left_point+1) - idx_couple->idx(left_point); //difference between the indices of two consecutive values of iP

        // update solution arrays
        if (sol_idx == -1){ // pre-computation not done
            for (int i = 0; i < num; i++){
                list_start_as_couple[i][idx] = tools::interp_1d_index_delta(par->grid_power, par->num_power, list_remain_couple[i], power_at_zero, left_point, delta, idx_couple->idx(0), 1, 0); 
            }
        }
        else{ // pre-computation done - get solution at sol_idx
            for (int i = 0; i < num; i++){
                list_start_as_couple[i][idx] = list_start_as_couple[i][idx_couple->idx(sol_idx)];
            }
        }
    } // end of update_to_indifference


    void check_participation_constraints(int* power_idx, double* power, double* Sw, double* Sm, index::index_couple_struct* idx_couple, double** list_start_as_couple, double** list_remain_couple, double* list_trans_to_single, int num, par_struct* par){

        // step 0: identify key indicators for each spouse
        // 0a: min and max surplus for each spouse
        double min_w = Sw[0];
        double max_w = Sw[par->num_power-1];
        double min_m = Sm[par->num_power-1];
        double max_m = Sm[0];

        // 0b: check if wife and husband have indifference points
        bool cross_w = (min_w < 0.0) && (max_w > 0.0);
        bool cross_m = (min_m < 0.0) && (max_m > 0.0);

        // 0b: check if wife and husband are always happy
        bool always_happy_w = (min_w > 0.0);
        bool always_happy_m = (min_m > 0.0);

        // 0c: check if wife and husband are never happy
        bool never_happy_w = (max_w < 0.0);
        bool never_happy_m = (max_m < 0.0);

        // step 1: check endpoints
        // 1a. check if all values are consistent with marriage
        if (always_happy_w && always_happy_m){
            for(int iP=0; iP<par->num_power; iP++){
                remain(iP, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par);
            }
        } //case 1a

        // 1b. check if all values are consistent with divorce
        else if (never_happy_w || never_happy_m){

            for (int iP=0; iP<par->num_power; iP++){
                divorce(iP, power_idx, power, idx_couple, list_start_as_couple, list_trans_to_single, num, par);
            }
        } //case 1b

        // 1c. check if husband is always happy, wife has indifference point
        else if (cross_w && always_happy_m){
            // find wife's indifference point
            int left_w = find_left(Sw, par->num_power);
            int Low_w = left_w+1;
            double power_at_zero_w = tools::interp_1d_index(Sw, par->num_power, par->grid_power, 0.0, left_w);

            // update case 1c
            for (int iP=0; iP<par->num_power; iP++){
                if (iP == 0){
                    update_to_indifference(iP, left_w, Low_w, power_at_zero_w, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par, -1);
                }
                else if (iP < Low_w){
                    update_to_indifference(iP, left_w, Low_w, power_at_zero_w, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par, 0);
                }
                else{
                    remain(iP, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par);
                } //if
            } //for
        } //case 1c

        // 1d: check if wife is always happy, husband has indifference point
        else if (cross_m && always_happy_w){
            //find husband's indifference point
            int left_m = find_left(Sm, par->num_power);
            int Low_m = left_m;
            double power_at_zero_m = tools::interp_1d_index(Sm, par->num_power, par->grid_power, 0.0, left_m);

            // update case 1d
            for (int iP=0; iP<par->num_power; iP++){
                if (iP<=Low_m){
                    remain(iP, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par);
                }
                else if (iP==Low_m+1){
                    update_to_indifference(iP, left_m, Low_m, power_at_zero_m, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par, -1);
                }
                else{
                    update_to_indifference(iP, left_m, Low_m, power_at_zero_m, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par, Low_m+1);
                } //if
            } //for
        } //case 1d

        // 1e: Both have indifference points
        else {
            //find indifference points
            int left_w = find_left(Sw, par->num_power);
            int Low_w = left_w+1;
            double power_at_zero_w = tools::interp_1d_index(Sw, par->num_power, par->grid_power, 0.0, left_w);

            int left_m = find_left(Sm, par->num_power);
            int Low_m = left_m;         
            double power_at_zero_m = tools::interp_1d_index(Sm, par->num_power, par->grid_power, 0.0, left_m);

            // update case 1e
            // no room for bargaining
            if (power_at_zero_w>power_at_zero_m) {
                for (int iP=0; iP<par->num_power; iP++){
                    divorce(iP, power_idx, power, idx_couple, list_start_as_couple, list_trans_to_single, num, par);
                }
            }
            //bargaining
            else {
                for (int iP=0; iP<par->num_power; iP++){
                    if (iP==0){ //update to woman's indifference point
                        update_to_indifference(iP, left_w, Low_w, power_at_zero_w, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par, -1);
                    }
                    else if (iP<Low_w){ //re-use pre-computed values
                        update_to_indifference(iP, left_w, Low_w, power_at_zero_w, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par, 0);
                    }
                    else if (iP>=Low_w && iP <= Low_m) { //no change between Low_w and Low_m
                        remain(iP, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par);
                    }
                    else if (iP == Low_m+1) { //update to man's indifference point
                        update_to_indifference(iP, left_m, Low_m, power_at_zero_m, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par, -1);
                    }
                    else { // re-use precomputed values
                        update_to_indifference(iP, left_m, Low_m, power_at_zero_m, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par, Low_m+1);
                    } //if (indifference points)
                }//for
            } //if (bargaining)
        } //case 1e
    } //end of check_participation_constraints


    void check_participation_constraints_verbose(int* power_idx, double* power, double* Sw, double* Sm, index::index_couple_struct* idx_couple, double** list_start_as_couple, double** list_remain_couple, double* list_trans_to_single, int num, par_struct* par, bool do_print=false){

        // step 0: identify key indicators for each spouse
        // 0a: min and max surplus for each spouse
        double min_w = Sw[0];
        double max_w = Sw[par->num_power-1];
        double min_m = Sm[par->num_power-1];
        double max_m = Sm[0];

        if (do_print) {
            logs::write("barg_log.txt", 1, "\n\nInitial checks");
            logs::write("barg_log.txt", 1, "\n   - min surplus for wife: %f", min_w);
            logs::write("barg_log.txt", 1, "\n   - max surplus for wife: %f", max_w);
            logs::write("barg_log.txt", 1, "\n   - min surplus for husband: %f", min_m);
            logs::write("barg_log.txt", 1, "\n   - max surplus for husband: %f", max_m);
        }

        // 0b: check if wife and husband have indifference points
        bool cross_w = (min_w < 0.0) && (max_w > 0.0);
        bool cross_m = (min_m < 0.0) && (max_m > 0.0);

        // 0b: check if wife and husband are always happy
        bool always_happy_w = (min_w > 0.0);
        bool always_happy_m = (min_m > 0.0);

        // 0c: check if wife and husband are never happy
        bool never_happy_w = (max_w < 0.0);
        bool never_happy_m = (max_m < 0.0);

        // step 1: check endpoints
        // 1a. check if all values are consistent with marriage
        if (always_happy_w && always_happy_m){
            if (do_print) {
                logs::write("barg_log.txt", 1, "\n\nCase 1a: all values are consistent with marriage");
                logs::write("barg_log.txt", 1, "\n   - remain married, no change in power");
            }
            for(int iP=0; iP<par->num_power; iP++){
                remain(iP, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par);
            }
        } //case 1a

        // 1b. check if all values are consistent with divorce
        else if (never_happy_w || never_happy_m){
            if (do_print) {
                logs::write("barg_log.txt", 1, "\nCase 1b: all values are consistent with divorce");
                logs::write("barg_log.txt", 1, "\n   - divorce");
            }

            for (int iP=0; iP<par->num_power; iP++){
                divorce(iP, power_idx, power, idx_couple, list_start_as_couple, list_trans_to_single, num, par);
            }
        } //case 1b

        // 1c. check if husband is always happy, wife has indifference point
        else if (cross_w && always_happy_m){
            // find wife's indifference point
            int left_w = find_left(Sw, par->num_power);
            int Low_w = left_w+1;
            double power_at_zero_w = tools::interp_1d_index(Sw, par->num_power, par->grid_power, 0.0, left_w);

            if (do_print) {
                logs::write("barg_log.txt", 1, "\nCase 1c: husband is always happy, wife has indifference point");
                logs::write("barg_log.txt", 1, "\n   - Wife's indifference point at index %d", Low_w);
                logs::write("barg_log.txt", 1, "\n   - Wife's power at indifference point: %f", power_at_zero_w);
            }

            // update case 1c
            for (int iP=0; iP<par->num_power; iP++){
                if (iP == 0){
                    update_to_indifference(iP, left_w, Low_w, power_at_zero_w, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par, -1);
                    if (do_print) {
                        logs::write("barg_log.txt", 1, "\n       - updating index %d to wife's indifference point %d", iP, Low_w);
                    }
                }
                else if (iP < Low_w){
                    update_to_indifference(iP, left_w, Low_w, power_at_zero_w, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par, 0);
                    if (do_print) {
                        logs::write("barg_log.txt", 1, "\n       - updating index %d to wife's indifference point %d", iP, Low_w);
                    }
                }
                else{
                    remain(iP, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par);
                    if (do_print) {
                        logs::write("barg_log.txt", 1, "\n       - index %d remains unchanged", iP);
                    }
                } //if
            } //for
        } //case 1c

        // 1d: check if wife is always happy, husband has indifference point
        else if (cross_m && always_happy_w){
            //find husband's indifference point
            int left_m = find_left(Sm, par->num_power);
            int Low_m = left_m;
            double power_at_zero_m = tools::interp_1d_index(Sm, par->num_power, par->grid_power, 0.0, left_m);

            if (do_print) {
                logs::write("barg_log.txt", 1, "\nCase 1d: wife is always happy, husband has indifference point");
                logs::write("barg_log.txt", 1, "\n   - Husband's indifference point at index %d", Low_m);
                logs::write("barg_log.txt", 1, "\n   - Husband's power at indifference point: %f", power_at_zero_m);
            }

            // update case 1d
            for (int iP=0; iP<par->num_power; iP++){
                if (iP<=Low_m){
                    remain(iP, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par);
                    if (do_print) {
                        logs::write("barg_log.txt", 1, "\n       - index %d remains unchanged", iP);
                    }
                }
                else if (iP==Low_m+1){
                    update_to_indifference(iP, left_m, Low_m, power_at_zero_m, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par, -1);
                    if (do_print) {
                        logs::write("barg_log.txt", 1, "\n       - updating index %d to husband's indifference point %d", iP, Low_m);
                    }
                }
                else{
                    update_to_indifference(iP, left_m, Low_m, power_at_zero_m, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par, Low_m+1);
                    if (do_print) {
                        logs::write("barg_log.txt", 1, "\n       - updating index %d to husband's indifference point %d", iP, Low_m);
                    }
                } //if
            } //for
        } //case 1d

        // 1e: Both have indifference points
        else {
            //find indifference points
            int left_w = find_left(Sw, par->num_power);
            int Low_w = left_w+1;
            double power_at_zero_w = tools::interp_1d_index(Sw, par->num_power, par->grid_power, 0.0, left_w);

            int left_m = find_left(Sm, par->num_power);
            int Low_m = left_m;         
            double power_at_zero_m = tools::interp_1d_index(Sm, par->num_power, par->grid_power, 0.0, left_m);

            if (do_print) {
                logs::write("barg_log.txt", 1, "\nCase 1e: Both have indifference points");
                logs::write("barg_log.txt", 1, "\n   - Wife's indifference point at index %d", Low_w);
                logs::write("barg_log.txt", 1, "\n   - Wife's power at indifference point: %f", power_at_zero_w);
                logs::write("barg_log.txt", 1, "\n   - Husband's indifference point at index %d", Low_m);
                logs::write("barg_log.txt", 1, "\n   - Husband's power at indifference point: %f", power_at_zero_m);
            }

            // update case 1e
            // no room for bargaining
            if (power_at_zero_w>power_at_zero_m) {
                for (int iP=0; iP<par->num_power; iP++){
                    divorce(iP, power_idx, power, idx_couple, list_start_as_couple, list_trans_to_single, num, par);
                    if (do_print){
                        logs::write("barg_log.txt", 1, "\n       - no value consistent with marriage -> divorce");
                    }
                }
            }
            //bargaining
            else {
                for (int iP=0; iP<par->num_power; iP++){
                    if (iP==0){ //update to woman's indifference point
                        update_to_indifference(iP, left_w, Low_w, power_at_zero_w, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par, -1);
                        if (do_print){
                            logs::write("barg_log.txt", 1, "\n       - updating index %d to wife's indifference point %d", iP, Low_w);
                        }
                    }
                    else if (iP<Low_w){ //re-use pre-computed values
                        update_to_indifference(iP, left_w, Low_w, power_at_zero_w, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par, 0);
                        if (do_print){
                            logs::write("barg_log.txt", 1, "\n       - updating index %d to wife's indifference point %d", iP, Low_w);
                        }
                    }
                    else if (iP>=Low_w && iP <= Low_m) { //no change between Low_w and Low_m
                        remain(iP, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par);
                        if (do_print){
                            logs::write("barg_log.txt", 1, "\n       - index %d remains unchanged", iP);
                        }
                    }
                    else if (iP == Low_m+1) { //update to man's indifference point
                        update_to_indifference(iP, left_m, Low_m, power_at_zero_m, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par, -1);
                        if (do_print){
                            logs::write("barg_log.txt", 1, "\n       - updating index %d to husband's indifference point %d", iP, Low_m);
                        }
                    }
                    else { // re-use precomputed values
                        update_to_indifference(iP, left_m, Low_m, power_at_zero_m, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, num, par, Low_m+1);
                        if (do_print){
                            logs::write("barg_log.txt", 1, "\n       - updating index %d to husband's indifference point %d", iP, Low_m);
                        }
                    } //if (indifference points)
                }//for
            } //if (bargaining)
        } //case 1e
    } //end of check_participation_constraints


    int initial_weight(double* Sw,double* Sm,par_struct* par){
    // determine the initial bargaining weight. On grid for simplicity.
    double min = HUGE_VAL;
    int idx_power = -1;
    double obj = 0.0;
    for(int iP=0;iP<par->num_power;iP++){

        // calculate objective function (to be minimized)
        if((Sw[iP]>0) & (Sm[iP]>0)){
            // obj = - sqrt(Sw[iP])*sqrt(Sm[iP]);
            obj = - (Sw[iP])*(Sm[iP]); // identical

            // update minimum
            if(obj<min){
                min = obj;
                idx_power = iP;
            }

        }

    }

    return idx_power;
    
}

//NASH BARGAINED INITIAL WEIGHT
    typedef struct {        
        par_struct *par;   
        sol_struct *sol;                                    
        double (*surplus_func)(double,index::state_couple_struct*,index::state_single_struct*,int,par_struct*,sol_struct*); //surplus func as function of power and state
        index::state_couple_struct *state_couple;                                       // state - tbc 
        index::state_single_struct *state_single_w;                                       // state - tbc   
        index::state_single_struct *state_single_m;                                       // state - tbc                              
    } nash_solver_struct;

    // compute negative nash surplus for given power + nash_struct
    double objfunc_nash_bargain(unsigned n, const double *x, double *grad, void* solver_data_in){
        // unpack
        nash_solver_struct* solver_data = (nash_solver_struct*) solver_data_in; 
        par_struct* par = solver_data->par;
        sol_struct* sol = solver_data->sol;
        double (*surplus_func)(double,index::state_couple_struct*, index::state_single_struct*,int,par_struct*,sol_struct*) = solver_data->surplus_func; //TODO: take continuous states as input as generically as possible

        // calculate individual supluses
        double Sw_x = surplus_func(x[0],solver_data->state_couple, solver_data->state_single_w, woman, par, sol);
        double Sm_x = surplus_func(x[0],solver_data->state_couple, solver_data->state_single_m, man, par, sol);

        return -(Sw_x*Sm_x); 
    }

    
    double nash_bargain(nash_solver_struct* nash_struct){
        // for a given couple idx, find the initial bargaining weight

        // unpack
        par_struct* par = nash_struct->par;
        sol_struct* sol = nash_struct->sol;

        // set up solver
        int const dim = 1;
        auto opt = nlopt_create(NLOPT_LN_BOBYQA, dim);
        nlopt_set_min_objective(opt, objfunc_nash_bargain, nash_struct);

        // set bounds
        double lb[dim], ub[dim];
        lb[0] = par->grid_power[0];
        ub[0] = par->grid_power[par->num_power-1];
        nlopt_set_lower_bounds(opt, lb);
        nlopt_set_upper_bounds(opt, ub);

        //optimize
        double minf = 0.0;
        double mu[dim];
        mu[0] = 0.5;
        nlopt_optimize(opt, mu, &minf);
        nlopt_destroy(opt);

        // check surplus is positive
        double (*surplus_func)(double,index::state_couple_struct*, index::state_single_struct*,int,par_struct*,sol_struct*) = nash_struct->surplus_func; //TODO: take continuous states as input as generically as possible
        index::state_couple_struct* state_couple = nash_struct->state_couple;
        index::state_single_struct* state_single_w = nash_struct->state_single_w;
        index::state_single_struct* state_single_m = nash_struct->state_single_m;
        double Sw = surplus_func(mu[0], state_couple, state_single_w, woman, par, sol);
        double Sm = surplus_func(mu[0], state_couple, state_single_m, man, par, sol);
        if ((Sw<0.0) |(Sm<0.0)){
            mu[0] = -1.0;
        }

        return mu[0];
    }

} // namespace bargaining