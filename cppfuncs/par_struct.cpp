typedef struct par_struct
{
 double R;
 double beta;
 double div_A_share;
 double inc_w;
 double inc_m;
 double rho_w;
 double rho_m;
 double alpha1_w;
 double alpha1_m;
 double alpha2_w;
 double alpha2_m;
 double phi_w;
 double phi_m;
 int T;
 int num_A;
 double max_A;
 int num_power;
 int num_love;
 double max_love;
 double sigma_love;
 int num_shock_love;
 int num_Ctot;
 double max_Ctot;
 int seed;
 int simT;
 int simN;
 bool do_cpp;
 int threads;
 double* grid_A;
 double* grid_A_single;
 double* grid_Aw;
 double* grid_Am;
 double* grid_power;
 double* grid_love;
 double* grid_shock_love;
 double* grid_weight_love;
 double* grid_Ctot;
} par_struct;

double get_double_par_struct(par_struct* x, char* name){

 if( strcmp(name,"R") == 0 ){ return x->R; }
 else if( strcmp(name,"beta") == 0 ){ return x->beta; }
 else if( strcmp(name,"div_A_share") == 0 ){ return x->div_A_share; }
 else if( strcmp(name,"inc_w") == 0 ){ return x->inc_w; }
 else if( strcmp(name,"inc_m") == 0 ){ return x->inc_m; }
 else if( strcmp(name,"rho_w") == 0 ){ return x->rho_w; }
 else if( strcmp(name,"rho_m") == 0 ){ return x->rho_m; }
 else if( strcmp(name,"alpha1_w") == 0 ){ return x->alpha1_w; }
 else if( strcmp(name,"alpha1_m") == 0 ){ return x->alpha1_m; }
 else if( strcmp(name,"alpha2_w") == 0 ){ return x->alpha2_w; }
 else if( strcmp(name,"alpha2_m") == 0 ){ return x->alpha2_m; }
 else if( strcmp(name,"phi_w") == 0 ){ return x->phi_w; }
 else if( strcmp(name,"phi_m") == 0 ){ return x->phi_m; }
 else if( strcmp(name,"max_A") == 0 ){ return x->max_A; }
 else if( strcmp(name,"max_love") == 0 ){ return x->max_love; }
 else if( strcmp(name,"sigma_love") == 0 ){ return x->sigma_love; }
 else if( strcmp(name,"max_Ctot") == 0 ){ return x->max_Ctot; }
 else {return NAN;}

}


int get_int_par_struct(par_struct* x, char* name){

 if( strcmp(name,"T") == 0 ){ return x->T; }
 else if( strcmp(name,"num_A") == 0 ){ return x->num_A; }
 else if( strcmp(name,"num_power") == 0 ){ return x->num_power; }
 else if( strcmp(name,"num_love") == 0 ){ return x->num_love; }
 else if( strcmp(name,"num_shock_love") == 0 ){ return x->num_shock_love; }
 else if( strcmp(name,"num_Ctot") == 0 ){ return x->num_Ctot; }
 else if( strcmp(name,"seed") == 0 ){ return x->seed; }
 else if( strcmp(name,"simT") == 0 ){ return x->simT; }
 else if( strcmp(name,"simN") == 0 ){ return x->simN; }
 else if( strcmp(name,"threads") == 0 ){ return x->threads; }
 else {return -9999;}

}


bool get_bool_par_struct(par_struct* x, char* name){

 if( strcmp(name,"do_cpp") == 0 ){ return x->do_cpp; }
 else {return false;}

}


double* get_double_p_par_struct(par_struct* x, char* name){

 if( strcmp(name,"grid_A") == 0 ){ return x->grid_A; }
 else if( strcmp(name,"grid_A_single") == 0 ){ return x->grid_A_single; }
 else if( strcmp(name,"grid_Aw") == 0 ){ return x->grid_Aw; }
 else if( strcmp(name,"grid_Am") == 0 ){ return x->grid_Am; }
 else if( strcmp(name,"grid_power") == 0 ){ return x->grid_power; }
 else if( strcmp(name,"grid_love") == 0 ){ return x->grid_love; }
 else if( strcmp(name,"grid_shock_love") == 0 ){ return x->grid_shock_love; }
 else if( strcmp(name,"grid_weight_love") == 0 ){ return x->grid_weight_love; }
 else if( strcmp(name,"grid_Ctot") == 0 ){ return x->grid_Ctot; }
 else {return NULL;}

}


