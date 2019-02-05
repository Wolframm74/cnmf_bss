#include "local_inc.hpp"

/*Expose*/
mat::fixed<N_static, K_static> complex_arg_m3_Pwrt_C_Vnk_nk_negative_numerator;
mat::fixed<N_static, K_static> complex_arg_m3_Pwrt_C_Vnk_nk_positive_denominator;

static mat::fixed<F_static, N_static> magnitude_model_FxN_mat_local;

static mat::fixed<N_static, K_static> complex_arg_m3_Pwrt_C_Vnk_nk_positive_den;

/*V_nk related matrices*/
static mat::fixed<F_static, L_static> V_nk_comp_outmat_FQ;
static mat::fixed<F_static, K_static> V_nk_comp_outmat_FG;
static mat::fixed<F_static, K_static> V_nk_comp_outmat_FK;
static mat::fixed<F_static, N_static> V_nk_comp_outmat_FN;
static mat::fixed<O_static, K_static> V_nk_comp_outmat_OK;

static mat::fixed<N_static, K_static> complex_arg_m3_Pwrt_C_Vng_nk_positive_den;

/*V_ng related matrices*/
/*static mat::fixed<F_static, L_static> V_ng_comp_outmat_FQ;
static mat::fixed<N_static, K_static> V_ng_comp_outmat_FG;*/
/*static mat::fixed<F_static, K_static> V_ng_comp_outmat_FK;*/
static mat::fixed<F_static, N_static> V_ng_comp_outmat_FN;
/*static mat::fixed<O_static, K_static> V_ng_comp_outmat_OK;*/

static void complex_argument_costfun_m3_V_update4_entry_fun2(int m_index_b, int m_index_a, cx_cube* W_fom_cx_p, mat* Z_ol_p, mat* Y_lk_p, cx_cube* expj_Phi_S_fkn_p, mat* T_fk_p, mat* V_nk_p, cx_cube* Xtilde_fnm_p, cx_cube* Xhat_fnm_p, cube* W_fom_p, cube* Xhat_low_fnm_p){

/*magnitude model computation*/
magnitude_model_FxN_mat_local=(*Xhat_low_fnm_p).slice(m_index_b)%(*Xhat_low_fnm_p).slice(m_index_a);

/*V_nk partial derivative related computation*/
V_nk_comp_outmat_FQ=(*W_fom_p).slice(m_index_a)*(*Z_ol_p);

V_nk_comp_outmat_FG=V_nk_comp_outmat_FQ*(*Y_lk_p);

V_nk_comp_outmat_FG=V_nk_comp_outmat_FG%(*T_fk_p);

V_nk_comp_outmat_FN=V_nk_comp_outmat_FG*trans(*V_nk_p);

V_nk_comp_outmat_OK=(*Z_ol_p)*(*Y_lk_p);

V_nk_comp_outmat_FK=(*W_fom_p).slice(m_index_b)*V_nk_comp_outmat_OK;

V_nk_comp_outmat_FK=V_nk_comp_outmat_FK%(*T_fk_p);

complex_arg_m3_Pwrt_C_Vnk_nk_positive_den=2*trans(V_nk_comp_outmat_FN)*(V_nk_comp_outmat_FK);

/*V_ng partial derivative related computation*/

V_ng_comp_outmat_FN=V_nk_comp_outmat_FK*trans(*V_nk_p);

complex_arg_m3_Pwrt_C_Vng_nk_positive_den=2*trans(V_ng_comp_outmat_FN)*V_nk_comp_outmat_FG;

complex_arg_m3_Pwrt_C_Vnk_nk_positive_denominator=complex_arg_m3_Pwrt_C_Vnk_nk_positive_den+complex_arg_m3_Pwrt_C_Vng_nk_positive_den;

complex_arg_m3_Pwrt_C_Vnk_nk_negative_numerator=complex_arg_m3_Pwrt_C_Vnk_nk_negative_numerator+complex_arg_m3_Pwrt_C_Vnk_nk_positive_denominator;

}

static cx_mat::fixed<F_static, O_static> kn_th_dummy_FO_mat_1;
static cx_mat::fixed<F_static, O_static> kn_th_dummy_FO_mat_2;

static cx_cube::fixed<F_static, O_static, L_static> accu_cube_FOL_1;
static cx_cube::fixed<F_static, O_static, L_static> accu_cube_FOL_2;

/*static cx_mat::fixed<F_static, N_static> G_ba_tilde_FxN_cx_mat_local;
static cx_mat::fixed<F_static, N_static> G_ba_hat_FxN_cx_mat_local;*/
static cx_mat::fixed<F_static, N_static> Error_FxN_cx_mat_local;

static void complex_argument_costfun_m3_V_update4_entry_fun1(int m_index_b, int m_index_a, cx_cube* W_fom_cx_p, mat* Z_ol_p, mat* Y_lk_p, cx_cube* expj_Phi_S_fkn_p, mat* T_fk_p, mat* V_nk_p, cx_cube* Xtilde_fnm_p, cx_cube* Xhat_fnm_p){

int k_index, n_index;
int l_index;

/*G_ba_tilde_FxN_cx_mat_local=(*Xtilde_fnm_p).slice(m_index_b)%conj((*Xtilde_fnm_p).slice(m_index_a));
G_ba_hat_FxN_cx_mat_local=(*Xhat_fnm_p).slice(m_index_b)%conj((*Xhat_fnm_p).slice(m_index_a));*/

Error_FxN_cx_mat_local=(*Xtilde_fnm_p).slice(m_index_b)%conj((*Xtilde_fnm_p).slice(m_index_a))-(*Xhat_fnm_p).slice(m_index_b)%conj((*Xhat_fnm_p).slice(m_index_a));

for (k_index=0; k_index<K_static; k_index++){

	for (n_index=0; n_index<N_static; n_index++){

		for (l_index=0; l_index<L_static; l_index++){

			kn_th_dummy_FO_mat_1=(*Y_lk_p)(l_index, k_index)*kron( (Error_FxN_cx_mat_local.col(n_index))%((*Xhat_fnm_p).slice(m_index_a).col(n_index))%((*T_fk_p).col(k_index))%(conj((*expj_Phi_S_fkn_p).slice(n_index).col(k_index))) , trans((*Z_ol_p).col(l_index)) );

			kn_th_dummy_FO_mat_2=(*Y_lk_p)(l_index, k_index)*kron( (Error_FxN_cx_mat_local.col(n_index))%(conj((*Xhat_fnm_p).slice(m_index_b).col(n_index)))%((*T_fk_p).col(k_index))%((*expj_Phi_S_fkn_p).slice(n_index).col(k_index)) , trans((*Z_ol_p).col(l_index)) );

			accu_cube_FOL_1.slice(l_index)=conj((*W_fom_cx_p).slice(m_index_b));

			accu_cube_FOL_1.slice(l_index)=accu_cube_FOL_1.slice(l_index)%kn_th_dummy_FO_mat_1;

			accu_cube_FOL_2.slice(l_index)=(*W_fom_cx_p).slice(m_index_a);

			accu_cube_FOL_2.slice(l_index)=accu_cube_FOL_2.slice(l_index)%kn_th_dummy_FO_mat_2;

		}

		complex_arg_m3_Pwrt_C_Vnk_nk_negative_numerator(n_index, k_index)=2*real(accu(accu_cube_FOL_1))+2*real(accu(accu_cube_FOL_2));

	}

}

}

void complex_argument_costfun_m3_V_update4_entry(arg_struct_t* argStruct_p){

int m_index;

complex_arg_m3_Pwrt_C_Vnk_nk_negative_numerator.zeros();
/*Pwrt_C_Vnk_nk_local_positive_denominator.zeros();*/

m_index=0; 

/*for (m_index=0; m_index<M_static; m_index++){*/

complex_argument_costfun_m3_V_update4_entry_fun1(1, 0, argStruct_p->W_fom_cx_p, argStruct_p->Z_ol_p, argStruct_p->Y_lk_p, argStruct_p->expj_Phi_S_fkn_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->Xtilde_fnm_p, argStruct_p->Xhat_fnm_p);

complex_argument_costfun_m3_V_update4_entry_fun2(1, 0, argStruct_p->W_fom_cx_p, argStruct_p->Z_ol_p, argStruct_p->Y_lk_p, argStruct_p->expj_Phi_S_fkn_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->Xtilde_fnm_p, argStruct_p->Xhat_fnm_p, argStruct_p->W_fom_p, argStruct_p->Xhat_low_fnm_p);

/*populate_G_H_fn_channel1(m_index, argStruct_p->Xhat_fnm_p);

populate_P_G1_wrt_Vnk_fnk();

populate_P_H1_wrt_Vnk_fnk(m_index, argStruct_p->Xhat_fnm_p);

populate_P_arg_Xh1_wrt_Vnk_fnk();*/

/*Pwrt_C_Vnk_nk=Pwrt_C_Vnk_nk+compute_Pwrt_C_Vnk_fnk(m_index, argStruct_p->Xhat_fnm_p, argStruct_p->Xtilde_fnm_p);*/

/*}*/

/*V_update2_Vnk(argStruct_p->V_nk_p);*/

}