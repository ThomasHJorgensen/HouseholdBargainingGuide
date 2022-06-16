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
}