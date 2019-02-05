#include "local_inc.hpp"

/*H_fk matrices*/
static cx_mat::fixed<F_static, K_static> dummy_mat_conj_H_fk_1_meq_a;
static cx_mat::fixed<F_static, K_static> dummy_mat_conj_H_fk_1_meq_b;

static cx_cube::fixed<F_static, N_static, K_static> Pwrt_C_phi_fnk_local;

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
cxarg_m2_populate_dummy_mat_conj_H_fk_1(dummy_mat_conj_H_fk_1_meq_a, m_index, W_fom_cx_p, Z_ol_p, Y_lk_p);

m_index=m_index_b;

/*bth channel*/
cxarg_m2_populate_dummy_mat_conj_H_fk_1(dummy_mat_conj_H_fk_1_meq_b, m_index, W_fom_cx_p, Z_ol_p, Y_lk_p);

for (k_index=0; k_index<K_static; k_index++){

Pwrt_C_phi_fnk_local.slice(k_index)=-conj(G_ba_FxN_cx_mat)%(Pwrt_C_phi_fnk_local.slice(k_index)+kron( (dummy_mat_conj_H_fk_1_meq_a.col(k_index) )%((*T_fk_p).col(k_index)) , trans((*V_nk_p).col(k_index)) )%conj(rotated_expj_Phi_S_FNK.slice(k_index)))%(*Xhat_fnm_p).slice(m_index_b);

Pwrt_C_phi_fnk_local.slice(k_index)=Pwrt_C_phi_fnk_local.slice(k_index)-G_ba_FxN_cx_mat%(Pwrt_C_phi_fnk_local.slice(k_index)+kron( (dummy_mat_conj_H_fk_1_meq_b.col(k_index) )%((*T_fk_p).col(k_index)) , trans((*V_nk_p).col(k_index)) )%conj(rotated_expj_Phi_S_FNK.slice(k_index)))%(*Xhat_fnm_p).slice(m_index_a);

}


}

static cube::fixed<F_static, K_static, N_static> absval_tensor_FKN;	
static cube::fixed<N_static, K_static, F_static> absval_tensor_NKF; 

/*static cx_cube::fixed<F_static, N_static, K_static> phase_diff_cx_FNK;*/
static double epsilon_Phi_S_value;

static void populate_phase_diff_FNK(cx_cube* expj_Phi_S_fkn_p, cx_cube* expj_Phi_S_nkf_p, mat* T_fk_p, mat* V_nk_p){

int f_index, n_index, k_index;

epsilon_Phi_S_value=0.0001;

for (f_index=0; f_index<F_static; f_index++){

	for (n_index=0; n_index<N_static; n_index++){

		for (k_index=0; k_index<K_static; k_index++){

			(*expj_Phi_S_fkn_p)(f_index, k_index, n_index)=(*expj_Phi_S_fkn_p)(f_index, k_index, n_index)+2*epsilon_Phi_S_value*Pwrt_C_phi_fnk_local(f_index, n_index, k_index);

			(*expj_Phi_S_nkf_p)(n_index, k_index, f_index)=(*expj_Phi_S_fkn_p)(f_index, k_index, n_index);

		}

	}

}

/*}

static void populate_phase_diff_FNK_asdfsad(cx_cube* expj_Phi_S_fkn_p, cx_cube* expj_Phi_S_nkf_p, mat* T_fk_p, mat* V_nk_p){

int f_index, n_index, k_index;*/

/*rotate*/
/*Rotating function defined in complex_argument_costfun_m1.hpp*/
rotate_expj_Phi_S_nkf_p_to_FNK(expj_Phi_S_nkf_p);
/*Can now safely access an up-to-date rotated_expj_Phi_S_FNK, from this point on*/

/*On the way out compute a rank 1 approximation at each k_index and update T_fk and V_nk at each k_index*/
mat U_svd;
vec s_svd;
mat V_svd;

for (k_index=0; k_index<K_static; k_index++){

	svd(U_svd, s_svd, V_svd, abs(rotated_expj_Phi_S_FNK.slice(k_index)) );

	(*T_fk_p).col(k_index)=abs(s_svd(0)*U_svd.col(0))%(*T_fk_p).col(k_index);

	(*V_nk_p).col(k_index)=abs(V_svd.col(0))%(*V_nk_p).col(k_index);

}

/*With the updated rank1 approximation updated into T_fk and V_nk, can now Divide out the absolute value to obtain unit magnitude in each FxNxK bin. */
absval_tensor_NKF=abs((*expj_Phi_S_nkf_p));
absval_tensor_FKN=abs((*expj_Phi_S_fkn_p));

(*expj_Phi_S_nkf_p).elem(find(absval_tensor_NKF==0)).fill(1);
(*expj_Phi_S_fkn_p).elem(find(absval_tensor_FKN==0)).fill(1);

(*expj_Phi_S_nkf_p)=(*expj_Phi_S_nkf_p)/abs((*expj_Phi_S_nkf_p));

//FKN update below
(*expj_Phi_S_fkn_p)=(*expj_Phi_S_fkn_p)/abs((*expj_Phi_S_fkn_p));

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

populate_phase_diff_FNK(argStruct_p->expj_Phi_S_fkn_p, argStruct_p->expj_Phi_S_nkf_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p);

/*Phi_S_update2_entry();*/

}
