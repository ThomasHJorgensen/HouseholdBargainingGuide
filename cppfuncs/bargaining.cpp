
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
   
    int find_left_point(double *S, parstruct *par){

        if (S[0] > S[-1]){
            S = - S;
        }

        return tools::binary_search(0, par-> num_power, S, 0.0);
    }
}