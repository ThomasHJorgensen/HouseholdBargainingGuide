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

    typedef struct{
            int t;
            int iL;
            int iA;
            par_struct *par; 
            int idx(int iP){
                    return index::index4(t,iP,iL,iA , par->T,par->num_power,par->num_love,par->num_A); 
            }
        
    } index_couple_struct;
}