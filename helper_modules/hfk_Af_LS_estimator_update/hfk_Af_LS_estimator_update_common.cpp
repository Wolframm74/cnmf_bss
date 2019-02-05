#include "local_inc.hpp"

mat::fixed<O_static, K_static> hfk_Af_LS_estimator_update_quantity_Z_ok;

cx_cube::fixed<F_static, K_static, M_static> h_fkm_shared_quantity;
static cx_mat::fixed<F_static*K_static, M_static> h_fk_m_local_quantity;
static cx_mat::fixed<M_static, F_static*K_static> h_m_fk_local_quantity;

cx_cube::fixed<M_static, F_static, K_static> hfk_Af_LS_estimator_update_quantity_h_mfk;
cx_cube::fixed<K_static, F_static, N_static> hfk_Af_LS_estimator_update_quantity_s_hat_kfn;
cx_cube::fixed<M_static, K_static, F_static> hfk_Af_LS_estimator_update_quantity_A_mkf;

void hfk_Af_LS_estimator_update_compute_Z_ok(mat* Z_ol_p, mat* Y_lk_p){

hfk_Af_LS_estimator_update_quantity_Z_ok=(*Z_ol_p)*(*Y_lk_p);

}

void hfk_Af_LS_estimator_update_compute_h_mfk(cx_cube* W_fom_cx_p){

int m_index;

for (m_index=0; m_index<M_static; m_index++){

h_fkm_shared_quantity.slice(m_index)=(*W_fom_cx_p).slice(m_index)*(hfk_Af_LS_estimator_update_quantity_Z_ok);

}

/*rotate it*/

/*copy*/
arrayops::copy(h_fk_m_local_quantity.memptr() , h_fkm_shared_quantity.memptr(), F_static*K_static*M_static);

/*strans()*/
h_m_fk_local_quantity=strans(h_fk_m_local_quantity);

/*copy into target*/
arrayops::copy(hfk_Af_LS_estimator_update_quantity_h_mfk.memptr(), h_m_fk_local_quantity.memptr(), M_static*F_static*K_static);

}

void  hfk_Af_LS_estimator_update_compute_s_hat(mat* T_fk_p, mat* V_nk_p, cx_cube* expj_Phi_S_nkf_p){

int f_index, n_index;

for (n_index=0; n_index<N_static; n_index++){

	for (f_index=0; f_index<F_static; f_index++){

	/*hfk_Af_LS_estimator_update_quantity_s_hat_kfn.slice(n_index).col(f_index)=trans((*T_fk_p).row(f_index))%trans((*V_nk_p).row(n_index))%(strans((*expj_Phi_S_nkf_p).slice(f_index).row(n_index)));*/
	hfk_Af_LS_estimator_update_quantity_s_hat_kfn.slice(n_index).col(f_index)=strans(((*T_fk_p).row(f_index))%((*V_nk_p).row(n_index))%((*expj_Phi_S_nkf_p).slice(f_index).row(n_index)));

	}

}


}

static cx_mat::fixed<M_static, K_static> A_quantity1_mk_local;
static cx_mat::fixed<K_static, K_static> A_quantity2_mk_local;

static cx_mat::fixed<F_static*N_static, M_static> Xtilde_fn_m_local;
static cx_mat::fixed<M_static, F_static*N_static> Xtilde_m_fn_local;
static cx_cube::fixed<M_static, F_static, N_static> Xtilde_mfn_local;

static void rotate_Xtilde(cx_cube* Xtilde_fnm_p){

/*copy*/
arrayops::copy(Xtilde_fn_m_local.memptr() , (*Xtilde_fnm_p).memptr(), F_static*N_static*M_static);

/*strans()*/
Xtilde_m_fn_local=strans(Xtilde_fn_m_local);

/*copy into target*/
arrayops::copy(Xtilde_mfn_local.memptr(), Xtilde_m_fn_local.memptr(), M_static*F_static*N_static);

}

void  hfk_Af_LS_estimator_update_compute_A(cx_cube* Xtilde_fnm_p){

int f_index, n_index;

rotate_Xtilde(Xtilde_fnm_p);

for (f_index=0; f_index<F_static; f_index++){

	A_quantity1_mk_local.zeros();
	A_quantity2_mk_local.zeros();

	for (n_index=0; n_index<N_static; n_index++){

	A_quantity1_mk_local=A_quantity1_mk_local+(Xtilde_mfn_local.slice(n_index).col(f_index))*trans(hfk_Af_LS_estimator_update_quantity_s_hat_kfn.slice(n_index).col(f_index));
	A_quantity2_mk_local=A_quantity2_mk_local+(hfk_Af_LS_estimator_update_quantity_s_hat_kfn.slice(n_index).col(f_index))*trans(hfk_Af_LS_estimator_update_quantity_s_hat_kfn.slice(n_index).col(f_index));

	}

	A_quantity1_mk_local=A_quantity1_mk_local/((cx_double)N_static);
	A_quantity2_mk_local=A_quantity2_mk_local/((cx_double)N_static);

	hfk_Af_LS_estimator_update_quantity_A_mkf.slice(f_index)=A_quantity1_mk_local*pinv(A_quantity2_mk_local);

}		

}

void compute_h_mfk_s_hat_Amkf_common(arg_struct_t* argStruct_p){

hfk_Af_LS_estimator_update_compute_Z_ok(argStruct_p->Z_ol_p, argStruct_p->Y_lk_p);

hfk_Af_LS_estimator_update_compute_h_mfk(argStruct_p->W_fom_cx_p);

 hfk_Af_LS_estimator_update_compute_s_hat(argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->expj_Phi_S_nkf_p);

 hfk_Af_LS_estimator_update_compute_A(argStruct_p->Xtilde_fnm_p);

}

/*Two extra rotating helper functions below*/

static cx_mat::fixed<F_static*O_static, M_static> W_fo_m_local;
static cx_mat::fixed<M_static, F_static*O_static> W_m_fo_local;
cx_cube::fixed<M_static, F_static, O_static> W_mfo_cx_shared;

void rotate_W_fom_cx(cx_cube* W_fom_cx_p){

/*copy*/
arrayops::copy(W_fo_m_local.memptr() , (*W_fom_cx_p).memptr(), F_static*O_static*M_static);

/*strans()*/
W_m_fo_local=strans(W_fo_m_local);

/*copy into target*/
arrayops::copy(W_mfo_cx_shared.memptr(), W_m_fo_local.memptr(), M_static*F_static*O_static);

}


/*static cx_mat::fixed<F_static*O_static, M_static> expj_Phi_W_fo_m_local;
static cx_mat::fixed<M_static, F_static*O_static> expj_Phi_W_m_fo_local;
cx_cube::fixed<M_static, F_static, O_static> expj_Phi_W_mfo_shared;

void rotate_expj_Phi_W_fom(cx_cube* expj_Phi_W_fom_p){

// copy
arrayops::copy(expj_Phi_W_fo_m_local.memptr() , (*Xtilde_fnm_p).memptr(), F_static*O_static*M_static);

// strans()
expj_Phi_W_m_fo_local=strans(expj_Phi_W_fo_m_local);

// copy into target
arrayops::copy(expj_Phi_W_mfo_shared.memptr(), expj_Phi_W_m_fo_local.memptr(), M_static*F_static*O_static);

}*/


static cx_mat::fixed<M_static*F_static, K_static> h_mf_k_local_quantity;
static cx_mat::fixed<K_static, M_static*F_static> h_k_mf_local_quantity;
static cx_cube::fixed<K_static, M_static, F_static> h_kmf_local_quantity;

/*output is KxMxF*/
static void rotate_h_mfk_into_output_kmf(void){

/*rotate it*/

/*copy*/
arrayops::copy(h_mf_k_local_quantity.memptr() , hfk_Af_LS_estimator_update_quantity_h_mfk.memptr(), F_static*K_static*M_static);

/*strans()*/
h_k_mf_local_quantity=strans(h_mf_k_local_quantity);

/*copy into target*/
arrayops::copy(h_kmf_local_quantity.memptr(), h_k_mf_local_quantity.memptr(), M_static*F_static*K_static);

}

/*cx_cube::fixed<M_static, F_static, K_static> hfk_Af_LS_estimator_update_quantity_h_mfk;
cx_cube::fixed<M_static, K_static, F_static> hfk_Af_LS_estimator_update_quantity_A_mkf;*/

double J_value_secondary_objective;

void compute_J(void){

int f_index;
J_value_secondary_objective=0;

rotate_h_mfk_into_output_kmf();

for (f_index=0; f_index<F_static; f_index++){

J_value_secondary_objective=J_value_secondary_objective+real(as_scalar(trace((strans(h_kmf_local_quantity.slice(f_index))-hfk_Af_LS_estimator_update_quantity_A_mkf.slice(f_index))*trans(strans(h_kmf_local_quantity.slice(f_index))-hfk_Af_LS_estimator_update_quantity_A_mkf.slice(f_index)))));

}

}	