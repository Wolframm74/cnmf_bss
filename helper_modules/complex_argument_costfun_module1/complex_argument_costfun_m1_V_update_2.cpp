#include "local_inc.hpp"

static cx_cube::fixed<F_static, N_static, K_static> P_Xh1_wrt_Vnk_fnk;

static void populate_P_Xh1_wrt_Vnk_fnk(int m_index, cx_cube* W_fom_cx_p, mat* Z_ol_p, mat* Y_lk_p, cx_cube* expj_Phi_S_fkn_p, mat* T_fk_p, mat* V_nk_p){

int k_index, n_index; 

populate_dummy_mat_H_fk_1(m_index, W_fom_cx_p, Z_ol_p, Y_lk_p);

for (k_index=0; k_index<K_static; k_index++){

	for (n_index=0; n_index<N_static; n_index++){

	P_Xh1_wrt_Vnk_fnk.slice(k_index).col(n_index)=(dummy_mat_H_fk_1.col(k_index)%(*T_fk_p).col(k_index))%((*expj_Phi_S_fkn_p).slice(n_index).col(k_index));

	}

}

}

/*

The following functions and their data used to reside here, but now moved to "common_module.cpp"

static void populate_G_H_fn_channel1(int m_index, cx_cube* Xhat_fnm_p);

static void populate_G_H_fn_channel2(int m_index, cx_cube* Xhat_fnm_p);

*/


/*Function local data*/
/*mat::fixed<F_static, N_static> dummy_mat_real_FN;*/

/*Function local data*/

static cube::fixed<F_static, N_static, K_static> P_G1_wrt_Vnk_fnk;

static void populate_P_G1_wrt_Vnk_fnk(void){

/*Take the imaginary part of Pwrt_Xh1_fnk*/
P_G1_wrt_Vnk_fnk=imag(P_Xh1_wrt_Vnk_fnk);


}

/*Function local data*/

static cube::fixed<F_static, N_static, K_static> P_H1_wrt_Vnk_fnk;

static cube::fixed<F_static, N_static, K_static> P_H1_wrt_Vnk_num_fnk;

static void populate_P_H1_wrt_Vnk_fnk(int m_index, cx_cube* Xhat_fnm_p){

int n_index, k_index;

dummy_mat_real_FN=real((*Xhat_fnm_p).slice(m_index));

for (k_index=0; k_index<K_static; k_index++){

	for (n_index=0; n_index<N_static; n_index++){

P_H1_wrt_Vnk_num_fnk.slice(k_index).col(n_index)=(dummy_mat_real_FN.col(n_index))%(real(P_Xh1_wrt_Vnk_fnk.slice(k_index).col(n_index)));

P_H1_wrt_Vnk_num_fnk.slice(k_index).col(n_index)=((dummy_mat_real_FN.col(n_index))%(imag(P_Xh1_wrt_Vnk_fnk.slice(k_index).col(n_index))))+(P_H1_wrt_Vnk_num_fnk.slice(k_index).col(n_index));

P_H1_wrt_Vnk_fnk.slice(k_index).col(n_index)=(P_H1_wrt_Vnk_num_fnk.slice(k_index).col(n_index))/(sqrt_denmat_1_fn.col(n_index));

	}	

}	

}

/*---------------------------------------------------------------------------------------------------*/

static cube::fixed<F_static, N_static, K_static> P_arg_Xh1_wrt_Vnk_fnk;

static void populate_P_arg_Xh1_wrt_Vnk_fnk(void){

int k_index;

for (k_index=0; k_index<K_static; k_index++){

P_arg_Xh1_wrt_Vnk_fnk.slice(k_index)=2*ones_mat_FxN/(square(G1_outmat_fn)+square(H1_outmat_fn));

P_arg_Xh1_wrt_Vnk_fnk.slice(k_index)=(P_arg_Xh1_wrt_Vnk_fnk.slice(k_index))%((P_G1_wrt_Vnk_fnk.slice(k_index)%H1_outmat_fn)-(G1_outmat_fn%P_H1_wrt_Vnk_fnk.slice(k_index)));

}

}


/*Function local data*/
static mat::fixed<N_static, K_static> Pwrt_C_Vnk_nk;
static mat::fixed<N_static, K_static> Pwrt_C_Vnk_nk_mth_component;

static mat& compute_Pwrt_C_Vnk_fnk(int m_index, cx_cube* Xhat_fnm_p, cx_cube* Xtilde_fnm_p){

int n_index, k_index;

mat& Pwrt_C_Vnk_nk_mth_component_ref=Pwrt_C_Vnk_nk_mth_component;

Pwrt_C_Vnk_nk_mth_component.zeros();

/*Elementise wise FxNxK product of the thing you computed with the previous function vs the difference between Pwrt_arg_Xh1_fnk and Pwrt_arg_Xh2_fnk */
dummy_mat_real_FN.zeros();

dummy_mat_real_FN=compute_arg_X_FN(Xhat_fnm_p, m_index);

/*dummy_mat_real_FN=dummy_mat_real_FN+compute_arg_X_FN(Xtilde_fnm_p, m_index_2);*/

/*dummy_mat_real_FN=dummy_mat_real_FN-compute_arg_X_FN(Xtilde_fnm_p, m_index);*/

/*dummy_mat_real_FN=dummy_mat_real_FN-compute_arg_X_FN(Xhat_fnm_p, m_index_2);*/

dummy_mat_real_FN=2*dummy_mat_real_FN;

for (n_index=0; n_index<N_static; n_index++){

	for (k_index=0; k_index<K_static; k_index++){

		Pwrt_C_Vnk_nk_mth_component(n_index, k_index)=accu( dummy_mat_real_FN.col(n_index) % (P_arg_Xh1_wrt_Vnk_fnk.slice(k_index).col(n_index)) );

}

}

return Pwrt_C_Vnk_nk_mth_component_ref;

}

static void V_update2_Vnk(mat* V_nk_p){

double epsilon_Vnk;

epsilon_Vnk=1000;

(*V_nk_p)=(*V_nk_p)-(epsilon_Vnk)*(Pwrt_C_Vnk_nk);

}

void complex_argument_costfun_m1_V_update2_entry(arg_struct_t* argStruct_p){

int m_index;

Pwrt_C_Vnk_nk.zeros();

m_index=0; 

/*for (m_index=0; m_index<M_static; m_index++){*/

populate_P_Xh1_wrt_Vnk_fnk(m_index, argStruct_p->W_fom_cx_p, argStruct_p->Z_ol_p, argStruct_p->Y_lk_p, argStruct_p->expj_Phi_S_fkn_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p);

populate_G_H_fn_channel1(m_index, argStruct_p->Xhat_fnm_p);

populate_P_G1_wrt_Vnk_fnk();

populate_P_H1_wrt_Vnk_fnk(m_index, argStruct_p->Xhat_fnm_p);

populate_P_arg_Xh1_wrt_Vnk_fnk();

Pwrt_C_Vnk_nk=Pwrt_C_Vnk_nk+compute_Pwrt_C_Vnk_fnk(m_index, argStruct_p->Xhat_fnm_p, argStruct_p->Xtilde_fnm_p);

/*}*/

V_update2_Vnk(argStruct_p->V_nk_p);

}