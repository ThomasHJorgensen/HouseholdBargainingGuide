// functions related to the institutional environment.
#ifndef MAIN
#define TOOLS
#include "myheader.cpp"
#endif


namespace tools {

double maxf(double* vec, int num){

    double max = vec[0];
    for(int i=1;i<num;i++){
        if(vec[i]>max){
            max = vec[i];
        }
    }

    return max;
}
double minf(double* vec, int num){

    double min = vec[0];
    for(int i=1;i<num;i++){
        if(vec[i]<min){
            min = vec[i];
        }
    }

    return min;
}

int binary_search(int imin, int Nx, double *x, double xi)
{
    int imid, half;

    // a. checks
    if(xi <= x[0]){
        return 0;
    } else if(xi >= x[Nx-2]) {
        return Nx-2;
    }

    // b. binary search
    while((half = Nx/2)){
        imid = imin + half;
        imin = (x[imid] <= xi) ? imid:imin;
        Nx  -= half;
    }

    return imin;

}

double interp_1d_index(double* grid1,int num1 ,double* value1,double xi1,int j1){
    /* 1d interpolation for one point
        
    Args:
        grid1 : 1d grid
        value : value array (2d)
        xi1 : input point
    Returns:
        yi : output
    */

    // a. left/right
    double nom_left = grid1[j1+1]-xi1;
    double nom_right = xi1-grid1[j1];

    // b. interpolation
    double denom = (grid1[j1+1]-grid1[j1]);
    double nom = 0.0;
    for (size_t k1 = 0; k1 < 2; k1++){
        double nom_1 = nom_left;
        if (k1==1){
            nom_1 = nom_right;
        }
        nom += nom_1*value1[j1+k1];
    }

    return nom/denom;

} // interp_1d

double interp_1d(double* grid1,int num1 ,double* value1,double xi1){
    /* 1d interpolation for one point
        
    Args:
        grid1 : 1d grid
        value : value array (2d)
        xi1 : input point
    Returns:
        yi : output
    */

    // a. search 
    int j1 = binary_search(0,num1,grid1,xi1);

    return interp_1d_index(grid1,num1 ,value1,xi1,j1);

} // interp_1d

void interp_1d_2out_index(double* grid1,int num1 ,double* value1,double* value2,double xi1 , double* out1,double* out2 ,int j1){

    // a. left/right
    double nom_1_left = grid1[j1+1]-xi1;
    double nom_1_right = xi1-grid1[j1];

    // b. interpolation
    double denom = (grid1[j1+1]-grid1[j1]);
    double nom1 = 0.0;
    double nom2 = 0.0;
    for (size_t k1 = 0; k1 < 2; k1++){
        double nom_1 = nom_1_left;
        if (k1==1){
            nom_1 = nom_1_right;
        }
        nom1 += nom_1*value1[j1+k1];
        nom2 += nom_1*value2[j1+k1];
    }

    out1[0] = nom1/denom;
    out2[0] = nom2/denom;
}


void interp_1d_2out(double* grid1,int num1 ,double* value1,double* value2,double xi1 , double* out1,double* out2 ){
    /* 1d interpolation for one point
        
    Args:
        grid1 : 1d grid
        value : value array (2d)
        xi1 : input point
    Returns:
        yi : output
    */

    // a. search in each dimension
    int j1 = binary_search(0,num1,grid1,xi1);

    // b. calculate the interpolants using the index
    interp_1d_2out_index(grid1,num1,value1,value2,xi1,out1,out2 ,j1);

} // interp_1d_2out


void interp_2d_2out_index(double* grid1,double* grid2,int num2, double* value1,double* value2,double xi1,double xi2,int j1,int j2, double* out1,double* out2){
    // b. pre-calculate nom parts and denominator
    double nom_1_left = grid1[j1+1]-xi1;
    double nom_1_right = xi1-grid1[j1];

    double nom_2_left = grid2[j2+1]-xi2;
    double nom_2_right = xi2-grid2[j2];

    double denom = (grid1[j1+1]-grid1[j1])*(grid2[j2+1]-grid2[j2]);

    // c. interpolation
    double nom1 = 0.0;
    double nom2 = 0.0;
    for (size_t k1 = 0; k1 < 2; k1++){
        double nom_1 = nom_1_left;
        if (k1==1){ nom_1 = nom_1_right;}

        for (size_t k2 = 0; k2 < 2; k2++){
            double nom_2 = nom_2_left;
            if (k2==1){ nom_2 = nom_2_right;}

            nom1 += nom_1*nom_2*value1[(j1+k1)*num2 + j2+k2];
            nom2 += nom_1*nom_2*value2[(j1+k1)*num2 + j2+k2];
        }
    }
    
    out1[0] = nom1/denom;
    out2[0] = nom2/denom;

}

void interp_2d_2out(double* grid1,double* grid2,int num1, int num2,double* value1,double* value2,double xi1,double xi2, double* out1,double* out2){    
  
    // a. search in each dimension
    int j1 = binary_search(0,num1,grid1,xi1);
    int j2 = binary_search(0,num2,grid2,xi2);

    interp_2d_2out_index(grid1,grid2,num2,value1,value2,xi1,xi2,j1,j2,out1,out2);
}


double interp_2d(double* grid1,double* grid2,int num1, int num2,double* value,double xi1,double xi2){
    
    // a. search in each dimension
    int j1 = binary_search(0,num1,grid1,xi1);
    int j2 = binary_search(0,num2,grid2,xi2);

    // b. left/right
    double nom_1_left = grid1[j1+1]-xi1;
    double nom_1_right = xi1-grid1[j1];

    double nom_2_left = grid2[j2+1]-xi2;
    double nom_2_right = xi2-grid2[j2];

    // c. interpolation
    double denom = (grid1[j1+1]-grid1[j1])*(grid2[j2+1]-grid2[j2]);

    double nom = 0.0;
    for (size_t k1 = 0; k1 < 2; k1++){
        double nom_1 = nom_1_left;
        if (k1==1){ nom_1 = nom_1_right;}

        for (size_t k2 = 0; k2 < 2; k2++){
            double nom_2 = nom_2_left;
            if (k2==1){ nom_2 = nom_2_right;}

            int idx = index::index2(j1+k1,j2+k2,num1,num2);
            nom += nom_1*nom_2*value[idx]; //value[(j1+k1)*num2 + j2+k2];
        }
    }

    return nom/denom;
}


double interp_2d_int(double* grid1,double* grid2,int num1, int num2,int* value,double xi1,double xi2){
    
    // a. search in each dimension
    int j1 = binary_search(0,num1,grid1,xi1);
    int j2 = binary_search(0,num2,grid2,xi2);

    // b. left/right
    double nom_1_left = grid1[j1+1]-xi1;
    double nom_1_right = xi1-grid1[j1];

    double nom_2_left = grid2[j2+1]-xi2;
    double nom_2_right = xi2-grid2[j2];

    // c. interpolation
    double denom = (grid1[j1+1]-grid1[j1])*(grid2[j2+1]-grid2[j2]);

    double nom = 0.0;
    for (size_t k1 = 0; k1 < 2; k1++){
        double nom_1 = nom_1_left;
        if (k1==1){ nom_1 = nom_1_right;}

        for (size_t k2 = 0; k2 < 2; k2++){
            double nom_2 = nom_2_left;
            if (k2==1){ nom_2 = nom_2_right;}

            int idx = index::index2(j1+k1,j2+k2,num1,num2);
            nom += nom_1*nom_2*(double)value[idx]; //value[(j1+k1)*num2 + j2+k2];
        }
    }

    return nom/denom;
}

double _interp_3d(double* grid1,double* grid2,double* grid3,int num1, int num2, int num3,double* value,double xi1,double xi2,double xi3){

    // a. search in each dimension
    int j1 = binary_search(0,num1,grid1,xi1);
    int j2 = binary_search(0,num2,grid2,xi2);
    int j3 = binary_search(0,num3,grid3,xi3);

    // a. left/right
    double nom_1_left = grid1[j1+1]-xi1;
    double nom_1_right = xi1-grid1[j1];

    double nom_2_left = grid2[j2+1]-xi2;
    double nom_2_right = xi2-grid2[j2];

    double nom_3_left = grid3[j3+1]-xi3;
    double nom_3_right = xi3-grid3[j3];

    // b. interpolation
    double denom = (grid1[j1+1]-grid1[j1])*(grid2[j2+1]-grid2[j2])*(grid3[j3+1]-grid3[j3]);
    double nom = 0.0;
    for (size_t k1 = 0; k1 < 2; k1++){
        double nom_1 = nom_1_left;
        if (k1==1){ nom_1 = nom_1_right;}

        for (size_t k2 = 0; k2 < 2; k2++){
            double nom_2 = nom_2_left;
            if (k2==1){ nom_2 = nom_2_right;}

            for (size_t k3 = 0; k3 < 2; k3++){
                double nom_3 = nom_3_left;
                if (k3==1){ nom_3 = nom_3_right;}            

                int idx = index::index3(j1+k1,j2+k2,j3+k3, num1,num2, num3);   
                nom += nom_1*nom_2*nom_3*value[idx];
            }
        }
    }

    return nom/denom;
}



double interp_3d(double* grid1,double* grid2,double* grid3,int num1, int num2, int num3,double* value,double xi1,double xi2,double xi3){

    // check if first dimension only has one element. If so use linear_interp_2d from consav package
    // This handles that child capital might not be kept track of
    if (num1 == 1){
        return interp_2d(grid2,grid3,num2,num3,value,xi2,xi3);
    } else {
        return _interp_3d(grid1,grid2,grid3,num1,num2,num3,value,xi1,xi2,xi3);
    }
}


double _interp_4d(double* grid1,double* grid2,double* grid3,double* grid4,int num1, int num2, int num3, int num4,double* value,double xi1,double xi2,double xi3,double xi4){
    // a. search in each dimension
    int j1 = binary_search(0,num1,grid1,xi1);
    int j2 = binary_search(0,num2,grid2,xi2);
    int j3 = binary_search(0,num3,grid3,xi3);
    int j4 = binary_search(0,num4,grid4,xi4);

    // b. left/right
    double nom_1_left = grid1[j1+1]-xi1;
    double nom_1_right = xi1-grid1[j1];

    double nom_2_left = grid2[j2+1]-xi2;
    double nom_2_right = xi2-grid2[j2];

    double nom_3_left = grid3[j3+1]-xi3;
    double nom_3_right = xi3-grid3[j3];

    double nom_4_left = grid4[j4+1]-xi4;
    double nom_4_right = xi4-grid4[j4];

    // c. interpolation
    double denom = (grid1[j1+1]-grid1[j1])*(grid2[j2+1]-grid2[j2])*(grid3[j3+1]-grid3[j3])*(grid4[j4+1]-grid4[j4]);
    double nom = 0.0;
    for (size_t k1 = 0; k1 < 2; k1++){
        double nom_1 = nom_1_left;
        if (k1==1){ nom_1 = nom_1_right;}

        for (size_t k2 = 0; k2 < 2; k2++){
            double nom_2 = nom_2_left;
            if (k2==1){ nom_2 = nom_2_right;}

            for (size_t k3 = 0; k3 < 2; k3++){
                double nom_3 = nom_3_left;
                if (k3==1){ nom_3 = nom_3_right;}    

                for (size_t k4 = 0; k4 < 2; k4++){
                    double nom_4 = nom_4_left;
                    if (k4==1){ nom_4 = nom_4_right;}         

                    int idx = index::index4(j1+k1,j2+k2,j3+k3,j4+k4, num1,num2, num3,num4);   
                    nom += nom_1*nom_2*nom_3*nom_4*value[idx];
                }
            }
        }
    }

    return nom/denom;
}

double interp_4d(double* grid1,double* grid2,double* grid3,double* grid4,int num1, int num2, int num3, int num4,double* value,double xi1,double xi2,double xi3,double xi4){

    // check if first dimension only has one element. If so use linear_interp_2d from consav package
    // This handles that child capital might not be kept track of
    if (num1 == 1){
        return interp_3d(grid2,grid3,grid4,num2,num3,num4,value,xi2,xi3,xi4);
    } else {
        return _interp_4d(grid1,grid2,grid3,grid4,num1,num2,num3,num4,value,xi1,xi2,xi3,xi4);
    }
}


} // tools