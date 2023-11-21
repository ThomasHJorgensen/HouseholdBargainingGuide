typedef struct sol_struct
{
 double* Vw_single;
 double* Vm_single;
 double* Cw_priv_single;
 double* Cm_priv_single;
 double* Cw_pub_single;
 double* Cm_pub_single;
 double* Cw_tot_single;
 double* Cm_tot_single;
 double* Vw_trans_single;
 double* Vm_trans_single;
 double* Cw_priv_trans_single;
 double* Cm_priv_trans_single;
 double* Cw_pub_trans_single;
 double* Cm_pub_trans_single;
 double* Cw_tot_trans_single;
 double* Cm_tot_trans_single;
 double* Vw_couple;
 double* Vm_couple;
 double* Cw_priv_couple;
 double* Cm_priv_couple;
 double* C_pub_couple;
 double* C_tot_couple;
 double* Vw_remain_couple;
 double* Vm_remain_couple;
 double* Cw_priv_remain_couple;
 double* Cm_priv_remain_couple;
 double* C_pub_remain_couple;
 double* C_tot_remain_couple;
 double* Sw;
 double* Sm;
 int* power_idx;
 double* power;
 double* savings_vec;
 double* Vw_plus_vec;
 double* Vm_plus_vec;
 double* marg_V_couple;
 double* marg_V_remain_couple;
 double* EmargU_pd;
 double* C_tot_pd;
 double* M_pd;
 double* marg_Vw_single;
 double* marg_Vm_single;
 double* EmargUw_single_pd;
 double* C_totw_single_pd;
 double* Mw_single_pd;
 double* EmargUm_single_pd;
 double* C_totm_single_pd;
 double* Mm_single_pd;
 double* pre_Ctot_Cw_priv;
 double* pre_Ctot_Cm_priv;
 double* pre_Ctot_C_pub;
 double* solution_time;
} sol_struct;

double* get_double_p_sol_struct(sol_struct* x, char* name){

 if( strcmp(name,"Vw_single") == 0 ){ return x->Vw_single; }
 else if( strcmp(name,"Vm_single") == 0 ){ return x->Vm_single; }
 else if( strcmp(name,"Cw_priv_single") == 0 ){ return x->Cw_priv_single; }
 else if( strcmp(name,"Cm_priv_single") == 0 ){ return x->Cm_priv_single; }
 else if( strcmp(name,"Cw_pub_single") == 0 ){ return x->Cw_pub_single; }
 else if( strcmp(name,"Cm_pub_single") == 0 ){ return x->Cm_pub_single; }
 else if( strcmp(name,"Cw_tot_single") == 0 ){ return x->Cw_tot_single; }
 else if( strcmp(name,"Cm_tot_single") == 0 ){ return x->Cm_tot_single; }
 else if( strcmp(name,"Vw_trans_single") == 0 ){ return x->Vw_trans_single; }
 else if( strcmp(name,"Vm_trans_single") == 0 ){ return x->Vm_trans_single; }
 else if( strcmp(name,"Cw_priv_trans_single") == 0 ){ return x->Cw_priv_trans_single; }
 else if( strcmp(name,"Cm_priv_trans_single") == 0 ){ return x->Cm_priv_trans_single; }
 else if( strcmp(name,"Cw_pub_trans_single") == 0 ){ return x->Cw_pub_trans_single; }
 else if( strcmp(name,"Cm_pub_trans_single") == 0 ){ return x->Cm_pub_trans_single; }
 else if( strcmp(name,"Cw_tot_trans_single") == 0 ){ return x->Cw_tot_trans_single; }
 else if( strcmp(name,"Cm_tot_trans_single") == 0 ){ return x->Cm_tot_trans_single; }
 else if( strcmp(name,"Vw_couple") == 0 ){ return x->Vw_couple; }
 else if( strcmp(name,"Vm_couple") == 0 ){ return x->Vm_couple; }
 else if( strcmp(name,"Cw_priv_couple") == 0 ){ return x->Cw_priv_couple; }
 else if( strcmp(name,"Cm_priv_couple") == 0 ){ return x->Cm_priv_couple; }
 else if( strcmp(name,"C_pub_couple") == 0 ){ return x->C_pub_couple; }
 else if( strcmp(name,"C_tot_couple") == 0 ){ return x->C_tot_couple; }
 else if( strcmp(name,"Vw_remain_couple") == 0 ){ return x->Vw_remain_couple; }
 else if( strcmp(name,"Vm_remain_couple") == 0 ){ return x->Vm_remain_couple; }
 else if( strcmp(name,"Cw_priv_remain_couple") == 0 ){ return x->Cw_priv_remain_couple; }
 else if( strcmp(name,"Cm_priv_remain_couple") == 0 ){ return x->Cm_priv_remain_couple; }
 else if( strcmp(name,"C_pub_remain_couple") == 0 ){ return x->C_pub_remain_couple; }
 else if( strcmp(name,"C_tot_remain_couple") == 0 ){ return x->C_tot_remain_couple; }
 else if( strcmp(name,"Sw") == 0 ){ return x->Sw; }
 else if( strcmp(name,"Sm") == 0 ){ return x->Sm; }
 else if( strcmp(name,"power") == 0 ){ return x->power; }
 else if( strcmp(name,"savings_vec") == 0 ){ return x->savings_vec; }
 else if( strcmp(name,"Vw_plus_vec") == 0 ){ return x->Vw_plus_vec; }
 else if( strcmp(name,"Vm_plus_vec") == 0 ){ return x->Vm_plus_vec; }
 else if( strcmp(name,"marg_V_couple") == 0 ){ return x->marg_V_couple; }
 else if( strcmp(name,"marg_V_remain_couple") == 0 ){ return x->marg_V_remain_couple; }
 else if( strcmp(name,"EmargU_pd") == 0 ){ return x->EmargU_pd; }
 else if( strcmp(name,"C_tot_pd") == 0 ){ return x->C_tot_pd; }
 else if( strcmp(name,"M_pd") == 0 ){ return x->M_pd; }
 else if( strcmp(name,"marg_Vw_single") == 0 ){ return x->marg_Vw_single; }
 else if( strcmp(name,"marg_Vm_single") == 0 ){ return x->marg_Vm_single; }
 else if( strcmp(name,"EmargUw_single_pd") == 0 ){ return x->EmargUw_single_pd; }
 else if( strcmp(name,"C_totw_single_pd") == 0 ){ return x->C_totw_single_pd; }
 else if( strcmp(name,"Mw_single_pd") == 0 ){ return x->Mw_single_pd; }
 else if( strcmp(name,"EmargUm_single_pd") == 0 ){ return x->EmargUm_single_pd; }
 else if( strcmp(name,"C_totm_single_pd") == 0 ){ return x->C_totm_single_pd; }
 else if( strcmp(name,"Mm_single_pd") == 0 ){ return x->Mm_single_pd; }
 else if( strcmp(name,"pre_Ctot_Cw_priv") == 0 ){ return x->pre_Ctot_Cw_priv; }
 else if( strcmp(name,"pre_Ctot_Cm_priv") == 0 ){ return x->pre_Ctot_Cm_priv; }
 else if( strcmp(name,"pre_Ctot_C_pub") == 0 ){ return x->pre_Ctot_C_pub; }
 else if( strcmp(name,"solution_time") == 0 ){ return x->solution_time; }
 else {return NULL;}

}


int* get_int_p_sol_struct(sol_struct* x, char* name){

 if( strcmp(name,"power_idx") == 0 ){ return x->power_idx; }
 else {return NULL;}

}


