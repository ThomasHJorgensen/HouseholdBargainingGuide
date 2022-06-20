typedef struct sim_struct
{
 double* Cw_priv;
 double* Cm_priv;
 double* Cw_pub;
 double* Cm_pub;
 double* Cw_tot;
 double* Cm_tot;
 double* C_tot;
 double* A;
 double* Aw;
 double* Am;
 double* couple;
 int* power_idx;
 double* power;
 double* love;
 double* draw_love;
 double* init_A;
 double* init_Aw;
 double* init_Am;
 bool* init_couple;
 int* init_power_idx;
 double* init_love;
} sim_struct;

double* get_double_p_sim_struct(sim_struct* x, char* name){

 if( strcmp(name,"Cw_priv") == 0 ){ return x->Cw_priv; }
 else if( strcmp(name,"Cm_priv") == 0 ){ return x->Cm_priv; }
 else if( strcmp(name,"Cw_pub") == 0 ){ return x->Cw_pub; }
 else if( strcmp(name,"Cm_pub") == 0 ){ return x->Cm_pub; }
 else if( strcmp(name,"Cw_tot") == 0 ){ return x->Cw_tot; }
 else if( strcmp(name,"Cm_tot") == 0 ){ return x->Cm_tot; }
 else if( strcmp(name,"C_tot") == 0 ){ return x->C_tot; }
 else if( strcmp(name,"A") == 0 ){ return x->A; }
 else if( strcmp(name,"Aw") == 0 ){ return x->Aw; }
 else if( strcmp(name,"Am") == 0 ){ return x->Am; }
 else if( strcmp(name,"couple") == 0 ){ return x->couple; }
 else if( strcmp(name,"power") == 0 ){ return x->power; }
 else if( strcmp(name,"love") == 0 ){ return x->love; }
 else if( strcmp(name,"draw_love") == 0 ){ return x->draw_love; }
 else if( strcmp(name,"init_A") == 0 ){ return x->init_A; }
 else if( strcmp(name,"init_Aw") == 0 ){ return x->init_Aw; }
 else if( strcmp(name,"init_Am") == 0 ){ return x->init_Am; }
 else if( strcmp(name,"init_love") == 0 ){ return x->init_love; }
 else {return NULL;}

}


int* get_int_p_sim_struct(sim_struct* x, char* name){

 if( strcmp(name,"power_idx") == 0 ){ return x->power_idx; }
 else if( strcmp(name,"init_power_idx") == 0 ){ return x->init_power_idx; }
 else {return NULL;}

}


bool* get_bool_p_sim_struct(sim_struct* x, char* name){

 if( strcmp(name,"init_couple") == 0 ){ return x->init_couple; }
 else {return NULL;}

}


