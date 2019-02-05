#include "local_inc.hpp"

static mat::fixed<N_static, K_static> Pwrt_C_Vnk_nk_local_negative_numerator;
static mat::fixed<N_static, K_static> Pwrt_C_Vnk_nk_local_positive_denominator;


static cx_cube::fixed<F_static, N_static, K_static> P_Xha_wrt_Vnk_fnk;
static cx_cube::fixed<F_static, N_static, K_static> P_conj_Xha_wrt_Vnk_fnk;
static cx_cube::fixed<F_static, N_static, K_static> P_Xhb_wrt_Vnk_fnk;
static cx_cube::fixed<F_static, N_static, K_static> P_conj_Xhb_wrt_Vnk_fnk;

/*H_fk matrices*/
static cx_mat::fixed<F_static, K_static> dummy_mat_H_fk_1_meq_a;
static cx_mat::fixed<F_static, K_static> dummy_mat_conj_H_fk_1_meq_a;
static cx_mat::fixed<F_static, K_static> dummy_mat_H_fk_1_meq_b;
static cx_mat::fixed<F_static, K_static> dummy_mat_conj_H_fk_1_meq_b;

static void populate_P_Xh1_wrt_Vnk_fnk_ba_and_conjs(int m_index_b, int m_index_a, cx_cube* W_fom_cx_p, mat* Z_ol_p, mat* Y_lk_p, cx_cube* expj_Phi_S_fkn_p, mat* T_fk_p, mat* V_nk_p, cx_cube* Xtilde_fnm_p, cx_cube* Xhat_fnm_p){

int k_index, n_index; 

int m_index;

m_index=m_index_a;

/*Compute G_ba */
G_ba_FxN_cx_mat=(*Xtilde_fnm_p).slice(m_index_b)%conj((*Xtilde_fnm_p).slice(m_index_a));

/*ath channel*/
cxarg_m2_populate_dummy_mat_H_fk_1(dummy_mat_H_fk_1_meq_a, m_index, W_fom_cx_p, Z_ol_p, Y_lk_p);
cxarg_m2_populate_dummy_mat_conj_H_fk_1(dummy_mat_conj_H_fk_1_meq_a, m_index, W_fom_cx_p, Z_ol_p, Y_lk_p);

m_index=m_index_b;

/*bth channel*/
cxarg_m2_populate_dummy_mat_H_fk_1(dummy_mat_H_fk_1_meq_b, m_index, W_fom_cx_p, Z_ol_p, Y_lk_p);
cxarg_m2_populate_dummy_mat_conj_H_fk_1(dummy_mat_conj_H_fk_1_meq_b, m_index, W_fom_cx_p, Z_ol_p, Y_lk_p);

for (k_index=0; k_index<K_static; k_index++){

	for (n_index=0; n_index<N_static; n_index++){

	P_Xha_wrt_Vnk_fnk.slice(k_index).col(n_index)=(dummy_mat_H_fk_1_meq_a.col(k_index)%(*T_fk_p).col(k_index))%((*expj_Phi_S_fkn_p).slice(n_index).col(k_index));

	P_conj_Xha_wrt_Vnk_fnk.slice(k_index).col(n_index)=(dummy_mat_conj_H_fk_1_meq_a.col(k_index)%(*T_fk_p).col(k_index))%conj((*expj_Phi_S_fkn_p).slice(n_index).col(k_index));

	P_Xhb_wrt_Vnk_fnk.slice(k_index).col(n_index)=(dummy_mat_H_fk_1_meq_b.col(k_index)%(*T_fk_p).col(k_index))%((*expj_Phi_S_fkn_p).slice(n_index).col(k_index));

	P_conj_Xhb_wrt_Vnk_fnk.slice(k_index).col(n_index)=(dummy_mat_conj_H_fk_1_meq_b.col(k_index)%(*T_fk_p).col(k_index))%conj((*expj_Phi_S_fkn_p).slice(n_index).col(k_index));

	Pwrt_C_Vnk_nk_local_negative_numerator(n_index, k_index)=-real(accu(G_ba_FxN_cx_mat.col(n_index)%((P_conj_Xhb_wrt_Vnk_fnk.slice(k_index).col(n_index)%conj((*Xhat_fnm_p).slice(m_index_a).col(n_index)))+(P_conj_Xha_wrt_Vnk_fnk.slice(k_index).col(n_index)%conj((*Xhat_fnm_p).slice(m_index_b).col(n_index))))));

	Pwrt_C_Vnk_nk_local_negative_numerator(n_index, k_index)=Pwrt_C_Vnk_nk_local_negative_numerator(n_index, k_index)-real(accu(conj(G_ba_FxN_cx_mat.col(n_index))%((P_Xhb_wrt_Vnk_fnk.slice(k_index).col(n_index)%((*Xhat_fnm_p).slice(m_index_a).col(n_index)))+(P_Xha_wrt_Vnk_fnk.slice(k_index).col(n_index)%((*Xhat_fnm_p).slice(m_index_b).col(n_index))))));

	Pwrt_C_Vnk_nk_local_positive_denominator(n_index, k_index)=Pwrt_C_Vnk_nk_local_positive_denominator(n_index, k_index)+2*real(accu(((*Xhat_fnm_p).slice(m_index_b).col(n_index))%(P_Xhb_wrt_Vnk_fnk.slice(k_index).col(n_index))%(square(abs((*Xhat_fnm_p).slice(m_index_a).col(n_index))))));

	Pwrt_C_Vnk_nk_local_positive_denominator(n_index, k_index)=Pwrt_C_Vnk_nk_local_positive_denominator(n_index, k_index)+2*real(accu(((*Xhat_fnm_p).slice(m_index_a).col(n_index))%(P_Xha_wrt_Vnk_fnk.slice(k_index).col(n_index))%(square(abs((*Xhat_fnm_p).slice(m_index_b).col(n_index))))));

	}

}

}


static void V_update2_Vnk(mat* V_nk_p){

double epsilon_Vnk;

epsilon_Vnk=0.0001;

/*(*V_nk_p)=(*V_nk_p)-(epsilon_Vnk)*(Pwrt_C_Vnk_nk_local);*/
(*V_nk_p)=(*V_nk_p)%(abs(Pwrt_C_Vnk_nk_local_negative_numerator)/Pwrt_C_Vnk_nk_local_positive_denominator);

}

void complex_argument_costfun_m2_V_update3_entry(arg_struct_t* argStruct_p){

int m_index;

Pwrt_C_Vnk_nk_local_negative_numerator.zeros();
Pwrt_C_Vnk_nk_local_positive_denominator.zeros();

m_index=0; 

/*for (m_index=0; m_index<M_static; m_index++){*/

populate_P_Xh1_wrt_Vnk_fnk_ba_and_conjs(1, 0, argStruct_p->W_fom_cx_p, argStruct_p->Z_ol_p, argStruct_p->Y_lk_p, argStruct_p->expj_Phi_S_fkn_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->Xtilde_fnm_p, argStruct_p->Xhat_fnm_p);

/*populate_G_H_fn_channel1(m_index, argStruct_p->Xhat_fnm_p);

populate_P_G1_wrt_Vnk_fnk();

populate_P_H1_wrt_Vnk_fnk(m_index, argStruct_p->Xhat_fnm_p);

populate_P_arg_Xh1_wrt_Vnk_fnk();*/

/*Pwrt_C_Vnk_nk=Pwrt_C_Vnk_nk+compute_Pwrt_C_Vnk_fnk(m_index, argStruct_p->Xhat_fnm_p, argStruct_p->Xtilde_fnm_p);*/

/*}*/

V_update2_Vnk(argStruct_p->V_nk_p);

}