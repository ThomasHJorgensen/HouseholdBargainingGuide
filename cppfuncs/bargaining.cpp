// functions for bargaining process.
#ifndef MAIN
#define BARGAINING
#include "myheader.cpp"
#endif

namespace bargaining {

    typedef struct{
        int t;
        int iL;
        int iA;
        par_struct *par; 

        int idx(int iP){
                return index::index4(t,iP,iL,iA , par->T,par->num_power,par->num_love,par->num_A); 
        }
    
    } index_couple_struct; //AMO: returns index for iP given t, iL and iA, like lambda function
   


    int find_left(double *S, par_struct *par){

        //flip is descending
        if (S[0] > S[-1]){
            for (int i = 0; i < par->num_power; i++){
                S[i] = -S[i];
            }
        }
        return tools::binary_search(0, par-> num_power, S, 0.0);
}

    // divorce
    void divorce(int iP, int *power_idx, double *power, int idx_single, index_couple_struct *idx_couple, double **list_start_as_couple, double **list_trans_to_single, int num, par_struct *par){
        int idx = idx_couple->idx(iP);
        power_idx[idx] = -1;
        power[idx] = -1.0;

        for (int i = 0; i < num; i++){
            list_start_as_couple[i][idx] = list_trans_to_single[i][idx_single];
        }
    }
    
    // remain
    void remain(int iP, int *power_idx, double *power, index_couple_struct *idx_couple, double **list_start_as_couple, double **list_remain_couple, int num, par_struct *par){
        int idx = idx_couple->idx(iP);
        power_idx[idx] = iP;
        power[idx] = par->grid_power[iP];



        for (int i = 0; i < num; i++){
            list_start_as_couple[i][idx] = list_remain_couple[i][iP];
        }
    }


    // update to indifference point
    void update_to_indifference(double *x, int iP, int left_point, double power_at_zero, int *power_idx, double *power, index_couple_struct *idx_couple, double **list_start_as_couple, double **list_remain_couple, int num, par_struct *par, int sol_idx = -1){
        int idx = idx_couple->idx(iP);
        power_idx[idx] = left_point;
        power[idx] = power_at_zero;


        // update solution arrays
        if (sol_idx == -1){ // pre-computation not done
            for (int i = 0; i < num; i++){
                list_start_as_couple[i][idx] = tools::interp_1d_index(par->grid_power, par->num_power, list_remain_couple[i], power_at_zero, left_point);
            }
        }
        else{ // pre-computation done - get solution at sol_idx
            for (int i = 0; i < num; i++){
                list_start_as_couple[i][idx] = list_start_as_couple[i][idx_couple->idx(sol_idx)];
            }
        }
    }

}