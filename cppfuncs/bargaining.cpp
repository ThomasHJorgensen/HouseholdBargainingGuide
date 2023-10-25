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


    void check_participation_constraints(int* power_idx, double* power, double* Sw, double* Sm, int* idx_single, index::index_couple_struct* idx_couple, double** list_start_as_couple, double** list_remain_couple, double* list_trans_to_single, int num, par_struct* par, bool do_print=false){

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
            double power_at_zero_w = tools::interp_1d_index(Sw, par->num_power, par->grid_power, 0.0, Low_w-1);

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
            double power_at_zero_m = tools::interp_1d_index(Sm, par->num_power, par->grid_power, 0.0, Low_m);

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
            double power_at_zero_w = tools::interp_1d_index(Sw, par->num_power, par->grid_power, 0.0, Low_w-1);

            int left_m = find_left(Sm, par->num_power);
            int Low_m = left_m;         
            double power_at_zero_m = tools::interp_1d_index(Sm, par->num_power, par->grid_power, 0.0, Low_m);

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

} // namespace bargaining