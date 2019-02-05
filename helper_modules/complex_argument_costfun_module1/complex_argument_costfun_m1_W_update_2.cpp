#include "local_inc.hpp"

/*get expj_Phi_S_nkf passed in*/

static rowvec::fixed<N_static> dummy_rowvec_real_f_1xN;
static rowvec::fixed<N_static> second_dummy_rowvec_real_f_1xN;

static cx_mat::fixed<N_static, K_static> dummy_mat_cx_NK;

static rowvec::fixed<K_static> dummy_rowvec_real_1xK;
static cx_rowvec::fixed<N_static> P_Xhm_wrt_Wfom_1xN;
static rowvec::fixed<N_static> P_Gm_wrt_Wfom_1xN;
static rowvec::fixed<N_static> P_Hm_wrt_Wfom_1xN;
static rowvec::fixed<N_static> P_Hm_wrt_Wfom_num_1xN;
static rowvec::fixed<N_static> P_arg_Xhm_wrt_Wfom_1xN;

static cube::fixed<F_static, O_static, M_static> Pwrt_C_Wfom_fom;

/*populate_G_H_f_1xN_channel_m local*/

static rowvec::fixed<N_static> Gm_outrow_f_1xN;
static rowvec::fixed<N_static> Hm_outrow_f_1xN;

static rowvec::fixed<N_static> sqrt_den_outrow_m_f_1xN;

static void W_update2_Wfom(cube* W_fom_p){

double epsilon_Wfom;

epsilon_Wfom=1000;

(*W_fom_p)=(*W_fom_p)-(epsilon_Wfom)*Pwrt_C_Wfom_fom;

}

/*static mat::fixed<F_static, N_static> dummy_mat_real_FN;*/

static void compute_Pwrt_C_Wfom_fom(cx_cube* expj_Phi_W_fom_p, mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, mat* Z_ol_p, cx_cube* Xhat_fnm_p, cx_cube* Xtilde_fnm_p, cx_cube* expj_Phi_S_nkf_p){

int m_index, f_index, o_index;

dummy_mat_real_FN.zeros();

/*for (m_index=0; m_index<M_static; m_index++){*/

	m_index=0; 

	dummy_mat_real_FN=2*compute_arg_X_FN(Xhat_fnm_p, m_index);

	/*dummy_mat_real_FN=dummy_mat_real_FN-2*compute_arg_X_FN(Xtilde_fnm_p, m_index);*/

	for (f_index=0; f_index<F_static; f_index++){

		/*populate G_H_channel_m_1xN*/

		Hm_outrow_f_1xN=real((*Xhat_fnm_p).slice(m_index).row(f_index));

		Gm_outrow_f_1xN=imag((*Xhat_fnm_p).slice(m_index).row(f_index));

		/*Compute  for the denominator function*/
		sqrt_den_outrow_m_f_1xN=sqrt(square(Gm_outrow_f_1xN)+square(Hm_outrow_f_1xN));
		/*May use the above result again later*/

		Hm_outrow_f_1xN=sqrt_den_outrow_m_f_1xN+Hm_outrow_f_1xN;

		second_dummy_rowvec_real_f_1xN=2*ones_row_1xN/(square(Hm_outrow_f_1xN)+square(Gm_outrow_f_1xN));

		/*end populate G_H_channel_m_1xN*/

		dummy_rowvec_real_f_1xN=real((*Xhat_fnm_p).slice(m_index).row(f_index));

		for (o_index=0; o_index<O_static; o_index++){

			dummy_mat_cx_NK=((*expj_Phi_W_fom_p)(f_index, o_index, m_index))*((*expj_Phi_S_nkf_p).slice(f_index));

			dummy_mat_cx_NK=dummy_mat_cx_NK%(*V_nk_p);

			dummy_rowvec_real_1xK=(*Z_ol_p).row(o_index)*(*Y_lk_p);

			dummy_rowvec_real_1xK=dummy_rowvec_real_1xK%((*T_fk_p).row(f_index));

			dummy_mat_cx_NK=dummy_mat_cx_NK%(kron(ones_col_Nx1, dummy_rowvec_real_1xK));

			P_Xhm_wrt_Wfom_1xN==trans(sum(dummy_mat_cx_NK, 1));

			P_Gm_wrt_Wfom_1xN=imag(P_Xhm_wrt_Wfom_1xN);	

			P_Hm_wrt_Wfom_num_1xN=(dummy_rowvec_real_f_1xN%real(P_Xhm_wrt_Wfom_1xN))+(Gm_outrow_f_1xN%imag(P_Xhm_wrt_Wfom_1xN));

			P_Hm_wrt_Wfom_1xN=P_Hm_wrt_Wfom_num_1xN/sqrt_den_outrow_m_f_1xN;	

			P_arg_Xhm_wrt_Wfom_1xN=second_dummy_rowvec_real_f_1xN%(P_Gm_wrt_Wfom_1xN%Hm_outrow_f_1xN+Gm_outrow_f_1xN%P_Hm_wrt_Wfom_1xN);

			Pwrt_C_Wfom_fom(f_index, o_index, m_index)=accu(dummy_mat_real_FN.row(f_index)%P_arg_Xhm_wrt_Wfom_1xN);

		}

	}

/*} */

}

void complex_argument_costfun_m1_W_update2_entry(arg_struct_t* argStruct_p){

compute_Pwrt_C_Wfom_fom(argStruct_p->expj_Phi_W_fom_p, argStruct_p->Y_lk_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->Z_ol_p, argStruct_p->Xhat_fnm_p, argStruct_p->Xtilde_fnm_p, argStruct_p->expj_Phi_S_nkf_p);

W_update2_Wfom(argStruct_p->W_fom_p);

}