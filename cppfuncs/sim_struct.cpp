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
 double* power;
 double* love;
 double* util;
 double* mean_lifetime_util;
 double* A_own;
 double* A_partner;
 double* draw_love;
 double* draw_meet;
 double* draw_uniform_partner_Aw;
 double* draw_uniform_partner_Am;
 int* draw_repartner_iL;
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
 else if( strcmp(name,"util") == 0 ){ return x->util; }
 else if( strcmp(name,"mean_lifetime_util") == 0 ){ return x->mean_lifetime_util; }
 else if( strcmp(name,"A_own") == 0 ){ return x->A_own; }
 else if( strcmp(name,"A_partner") == 0 ){ return x->A_partner; }
 else if( strcmp(name,"draw_love") == 0 ){ return x->draw_love; }
 else if( strcmp(name,"draw_meet") == 0 ){ return x->draw_meet; }
 else if( strcmp(name,"draw_uniform_partner_Aw") == 0 ){ return x->draw_uniform_partner_Aw; }
 else if( strcmp(name,"draw_uniform_partner_Am") == 0 ){ return x->draw_uniform_partner_Am; }
 else if( strcmp(name,"init_A") == 0 ){ return x->init_A; }
 else if( strcmp(name,"init_Aw") == 0 ){ return x->init_Aw; }
 else if( strcmp(name,"init_Am") == 0 ){ return x->init_Am; }
 else if( strcmp(name,"init_love") == 0 ){ return x->init_love; }
 else {return NULL;}

}


int* get_int_p_sim_struct(sim_struct* x, char* name){

 if( strcmp(name,"draw_repartner_iL") == 0 ){ return x->draw_repartner_iL; }
 else if( strcmp(name,"init_power_idx") == 0 ){ return x->init_power_idx; }
 else {return NULL;}

}


bool* get_bool_p_sim_struct(sim_struct* x, char* name){

 if( strcmp(name,"init_couple") == 0 ){ return x->init_couple; }
 else {return NULL;}

}


