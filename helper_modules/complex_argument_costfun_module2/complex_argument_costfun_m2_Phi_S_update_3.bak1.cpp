#include "local_inc.hpp"

/*Dummy mats*/
static cx_mat::fixed<F_static, N_static> dummy_mat_cx_FN_b;
static cx_mat::fixed<F_static, N_static> dummy_mat_cx_FN_conj_b;
static cx_mat::fixed<F_static, N_static> dummy_mat_cx_FN_a;
static cx_mat::fixed<F_static, N_static> dummy_mat_cx_FN_conj_a;

static cx_mat::fixed<F_static, N_static> dummy_mat_cx_FN_init_1;
static cx_mat::fixed<F_static, N_static> dummy_mat_cx_FN_init_conj_1;
/*mat::fixed<F_static, N_static> dummy_mat_real_FN;*/

/*m=1st channel*/

/*static cx_cube::fixed<F_static, N_static, K_static> Pwrt_Xh1_fnk;

static cx_cube::fixed<F_static, N_static, K_static> P_Xha_wrt_Phi_S_fnk;
static cx_cube::fixed<F_static, N_static, K_static> P_conj_Xha_wrt_Phi_S_fnk;
static cx_cube::fixed<F_static, N_static, K_static> P_Xhb_wrt_Phi_S_fnk;
static cx_cube::fixed<F_static, N_static, K_static> P_conj_Xhb_wrt_Phi_S_fnk;*/

/*H_fk matrices*/
static cx_mat::fixed<F_static, K_static> dummy_mat_H_fk_1_meq_a;
static cx_mat::fixed<F_static, K_static> dummy_mat_conj_H_fk_1_meq_a;
static cx_mat::fixed<F_static, K_static> dummy_mat_H_fk_1_meq_b;
static cx_mat::fixed<F_static, K_static> dummy_mat_conj_H_fk_1_meq_b;

static cube::fixed<F_static, N_static, K_static> Pwrt_C_phi_fnk_local;

static void populate_Pwrt_Xh1_fnk(int m_index_b, int m_index_a, cx_cube* W_fom_cx_p, mat* Z_ol_p, mat* Y_lk_p, cx_cube* expj_Phi_S_fkn_p, mat* T_fk_p, mat* V_nk_p, cx_cube* Xtilde_fnm_p, cx_cube* Xhat_fnm_p){
/*To answer the question why would you not just hardcore m_index to 1 for this function and 2 for the other: Because sometimes you might need to use m_indices not equal to just 1 and 2.
In the case that  you have >=2 microphones. 
*/

int k_index, n_index, m_index; 
cx_double fill_value, fill_value_conj;

fill_value.imag()=1; 

fill_value_conj.imag()=-1; 

/*Compute G_ba */
G_ba_FxN_cx_mat=(*Xtilde_fnm_p).slice(m_index_b)%conj((*Xtilde_fnm_p).slice(m_index_a));

m_index=m_index_a;

/*ath channel*/
cxarg_m2_populate_dummy_mat_H_fk_1(dummy_mat_H_fk_1_meq_a, m_index, W_fom_cx_p, Z_ol_p, Y_lk_p);
cxarg_m2_populate_dummy_mat_conj_H_fk_1(dummy_mat_conj_H_fk_1_meq_a, m_index, W_fom_cx_p, Z_ol_p, Y_lk_p);

m_index=m_index_b;

/*bth channel*/
cxarg_m2_populate_dummy_mat_H_fk_1(dummy_mat_H_fk_1_meq_b, m_index, W_fom_cx_p, Z_ol_p, Y_lk_p);
cxarg_m2_populate_dummy_mat_conj_H_fk_1(dummy_mat_conj_H_fk_1_meq_b, m_index, W_fom_cx_p, Z_ol_p, Y_lk_p);


for (k_index=0; k_index<K_static; k_index++){

	/*dummy_mat_cx_FN.fill(fill_value);*/

	dummy_mat_cx_FN_init_1.fill(fill_value);

	dummy_mat_cx_FN_init_conj_1.fill(fill_value_conj);

	dummy_mat_cx_FN_init_1=dummy_mat_cx_FN_init_1%((*T_fk_p).col(k_index)*trans((*V_nk_p).col(k_index)))%(rotated_expj_Phi_S_FNK.slice(k_index));

	dummy_mat_cx_FN_init_conj_1=dummy_mat_cx_FN_init_conj_1%((*T_fk_p).col(k_index)*trans((*V_nk_p).col(k_index)))%conj(rotated_expj_Phi_S_FNK.slice(k_index));

	/*ath channel*/
	dummy_mat_cx_FN_a=dummy_mat_cx_FN_init_1%(kron(dummy_mat_H_fk_1_meq_a.col(k_index), ones_row_1xN));

/*	dummy_mat_cx_FN_a=dummy_mat_cx_FN_a%((*T_fk_p).col(k_index)*trans((*V_nk_p).col(k_index)));*/

	dummy_mat_cx_FN_conj_a=dummy_mat_cx_FN_init_conj_1%(kron(dummy_mat_conj_H_fk_1_meq_a.col(k_index), ones_row_1xN));

/*	dummy_mat_cx_FN_conj_a=dummy_mat_cx_FN_conj_a%((*T_fk_p).col(k_index)*trans((*V_nk_p).col(k_index)));*/

	/*bth channel*/
	dummy_mat_cx_FN_b=dummy_mat_cx_FN_init_1%(kron(dummy_mat_H_fk_1_meq_b.col(k_index), ones_row_1xN));

/*	dummy_mat_cx_FN_b=dummy_mat_cx_FN_b%((*T_fk_p).col(k_index)*trans((*V_nk_p).col(k_index)));*/

	dummy_mat_cx_FN_conj_b=dummy_mat_cx_FN_init_conj_1%(kron(dummy_mat_conj_H_fk_1_meq_b.col(k_index), ones_row_1xN));

/*	dummy_mat_cx_FN_conj_b=dummy_mat_cx_FN_conj_b%((*T_fk_p).col(k_index)*trans((*V_nk_p).col(k_index)));*/

	/*Pwrt_Xh1_fnk.slice(k_index)=dummy_mat_cx_FN;*/

//	P_Xha_wrt_Phi_S_fnk.slice(k_index)=dummy_mat_cx_FN;

//	P_conj_Xha_wrt_Phi_S_fnk.slice(k_index)=dummy_mat_cx_FN_conj;

//	P_Xhb_wrt_Phi_S_fnk.slice(k_index)=dummy_mat_cx_FN;

//	P_conj_Xhb_wrt_Phi_S_fnk.slice(k_index)=dummy_mat_cx_FN_conj;

	/*for (n_index=0; n_index<N_static; n_index++){*/

	Pwrt_C_phi_fnk_local.slice(k_index)=-real(G_ba_FxN_cx_mat%((dummy_mat_cx_FN_conj_b%conj((*Xhat_fnm_p).slice(m_index_a)))+(dummy_mat_cx_FN_conj_a%conj((*Xhat_fnm_p).slice(m_index_b)))));

	Pwrt_C_phi_fnk_local.slice(k_index)=Pwrt_C_phi_fnk_local.slice(k_index)-real(conj(G_ba_FxN_cx_mat)%((dummy_mat_cx_FN_b%((*Xhat_fnm_p).slice(m_index_a)))+(dummy_mat_cx_FN_a%((*Xhat_fnm_p).slice(m_index_b)))));

	Pwrt_C_phi_fnk_local.slice(k_index)=Pwrt_C_phi_fnk_local.slice(k_index)+2*real(((*Xhat_fnm_p).slice(m_index_b))%(dummy_mat_cx_FN_b)%(square(abs((*Xhat_fnm_p).slice(m_index_a)))));

	Pwrt_C_phi_fnk_local.slice(k_index)=Pwrt_C_phi_fnk_local.slice(k_index)+2*real(((*Xhat_fnm_p).slice(m_index_a))%(dummy_mat_cx_FN_a)%(square(abs((*Xhat_fnm_p).slice(m_index_b)))));

	/*}*/

}


}



static cx_cube::fixed<F_static, N_static, K_static> phase_diff_cx_FNK;
static double epsilon_Phi_S_value;

static void populate_phase_diff_FNK(cx_cube* expj_Phi_S_fkn_p, cx_cube* expj_Phi_S_nkf_p){

int f_index, n_index, k_index;

epsilon_Phi_S_value=1;

phase_diff_cx_FNK.set_real(cos(-epsilon_Phi_S_value*Pwrt_C_phi_fnk_local));

phase_diff_cx_FNK.set_imag(sin(-epsilon_Phi_S_value*Pwrt_C_phi_fnk_local));

for (f_index=0; f_index<F_static; f_index++){

	for (n_index=0; n_index<N_static; n_index++){

		for (k_index=0; k_index<K_static; k_index++){

		((*expj_Phi_S_fkn_p)(f_index, k_index, n_index))=((*expj_Phi_S_fkn_p)(f_index, k_index, n_index))*phase_diff_cx_FNK(f_index, n_index, k_index);

		((*expj_Phi_S_nkf_p)(n_index, k_index, f_index))=((*expj_Phi_S_nkf_p)(n_index, k_index, f_index))*phase_diff_cx_FNK(f_index, n_index, k_index);

		}

	}		

}



}


void complex_argument_costfun_m2_Phi_S_update3_entry(arg_struct_t* argStruct_p){

int m_index;

/*for (m_index=0; m_index<M_static; m_index++){*/

m_index=0; 

rotate_expj_Phi_S_nkf_p_to_FNK(argStruct_p->expj_Phi_S_nkf_p);

populate_Pwrt_Xh1_fnk(1, 0, argStruct_p->W_fom_cx_p, argStruct_p->Z_ol_p, argStruct_p->Y_lk_p, argStruct_p->expj_Phi_S_fkn_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->Xtilde_fnm_p, argStruct_p->Xhat_fnm_p);

/*populate_G_H_fn_channel1_local(m_index, argStruct_p->Xhat_fnm_p);*/

/*populate_Pwrt_G1_fnk();*/

/*populate_Pwrt_H1_fnk(m_index, argStruct_p->Xhat_fnm_p);*/

/*populate_Pwrt_arg_Xh1_fnk();*/

/*compute_Pwrt_C_phi_fnk(m_index, argStruct_p->Xhat_fnm_p, argStruct_p->Xtilde_fnm_p);*/

/*}*/

populate_phase_diff_FNK(argStruct_p->expj_Phi_S_fkn_p, argStruct_p->expj_Phi_S_nkf_p);

/*Phi_S_update2_entry();*/

}
