// functions for calculating linear indices.
#ifndef MAIN
#define INDEX
#include "myheader.cpp"
#endif

namespace index {
    int index2(int i1,int i2,int N1,int N2){
        return i2 + i1*N2;
    }
    int index3(int i1,int i2,int i3,int N1,int N2, int N3){
        return i3 + i2*N3 + i1*N2*N3;
    }
    int index4(int i1,int i2,int i3,int i4,int N1,int N2, int N3, int N4){
        return i4 + (i3 + (i2 + i1*N2)*N3)*N4;
    }

    int couple(int t,int iP,int iL,int iA,par_struct* par){
        return index4(t,iP,iL,iA , par->T,par->num_power,par->num_love,par->num_A); 
    }
    int couple_pd(int t,int iP,int iL,int iA_pd,par_struct* par){
        return index4(t,iP,iL,iA_pd , par->T,par->num_power,par->num_love,par->num_A_pd); 
    }
    int single_to_couple(int t,int iL,int iA,par_struct* par){
        return index3(t,iL,iA , par->T,par->num_love,par->num_A); 
    }
    int single(int t,int iA,par_struct* par){
        return index2(t,iA , par->T,par->num_A); 
    }

    typedef struct{
            int t;
            int iL;
            int iA;
            par_struct *par; 
            int idx(int iP){
                    return index::couple(t,iP,iL,iA , par); 
            }
        
    } index_couple_struct;

    typedef struct{
            // state levels
            int t;
            double love;
            double A;
            double power;

            // indices
            int iL;
            int iA;

            // model content
            par_struct *par;
            sol_struct *sol;
    } state_couple_struct;

    typedef struct{
            // state levels
            int t;
            double A;

            // indices
            int iA;

            // model content
            par_struct *par;
            sol_struct *sol;
    } state_single_struct;
}