#include "local_inc.hpp"

cx_mat::fixed<F_static, N_static> G_ba_FxN_cx_mat; 

static cx_mat::fixed<N_static, K_static> dummy_mat_cx_NK;
static cx_mat::fixed<N_static, K_static> dummy_mat_cx_NK_conj;

static rowvec::fixed<K_static> dummy_rowvec_real_1xK;
static cx_rowvec::fixed<N_static> P_Xhm_wrt_Wfom_1xN;
static cx_rowvec::fixed<N_static> P_conj_Xhm_wrt_Wfom_1xN;

static cube::fixed<F_static, O_static, 2> Pwrt_C_Wfom_fom;

static colvec::fixed<2> m_index_array;

static void W_update3_Wfom(cube* W_fom_p){

double epsilon_Wfom;

epsilon_Wfom=0.1;

int m_iter, m_index; 

for (m_iter=0; m_iter<2; m_iter++){

m_index=m_index_array(m_iter);

(*W_fom_p).slice(m_index)=(*W_fom_p).slice(m_index)-(epsilon_Wfom)*Pwrt_C_Wfom_fom.slice(m_index);

}

}

/*static mat::fixed<F_static, N_static> dummy_mat_real_FN;*/



static void compute_Pwrt_C_Wfom_fom(int m_index_b, int m_index_a, cx_cube* expj_Phi_W_fom_p, mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, mat* Z_ol_p, cx_cube* Xhat_fnm_p, cx_cube* Xtilde_fnm_p, cx_cube* expj_Phi_S_nkf_p){

int m_index, f_index, o_index, m_index_other;
cx_double fill_value, fill_value_conj;

fill_value.imag()=1; 

fill_value_conj.imag()=-1; 

m_index_array(0)=m_index_b;
m_index_array(1)=m_index_a;

G_ba_FxN_cx_mat=(*Xtilde_fnm_p).slice(m_index_b)%conj((*Xtilde_fnm_p).slice(m_index_a));

int m_iter;


	for (f_index=0; f_index<F_static; f_index++){

		for (o_index=0; o_index<O_static; o_index++){

			for (m_iter=0; m_iter<2; m_iter++){

			m_index=m_index_array(m_iter);

			dummy_mat_cx_NK.fill(fill_value);

			dummy_mat_cx_NK_conj.fill(fill_value_conj);

			dummy_mat_cx_NK=((*expj_Phi_W_fom_p)(f_index, o_index, m_index))*((*expj_Phi_S_nkf_p).slice(f_index));

			dummy_mat_cx_NK_conj=conj((*expj_Phi_W_fom_p)(f_index, o_index, m_index))*conj((*expj_Phi_S_nkf_p).slice(f_index));

			dummy_mat_cx_NK=dummy_mat_cx_NK%(*V_nk_p);

			dummy_mat_cx_NK_conj=dummy_mat_cx_NK_conj%(*V_nk_p);

			dummy_rowvec_real_1xK=(*Z_ol_p).row(o_index)*(*Y_lk_p);

			dummy_rowvec_real_1xK=dummy_rowvec_real_1xK%((*T_fk_p).row(f_index));

			dummy_mat_cx_NK=dummy_mat_cx_NK%(kron(ones_col_Nx1, dummy_rowvec_real_1xK));

			dummy_mat_cx_NK_conj=dummy_mat_cx_NK_conj%(kron(ones_col_Nx1, dummy_rowvec_real_1xK));

			P_Xhm_wrt_Wfom_1xN==trans(sum(dummy_mat_cx_NK, 1));

			P_conj_Xhm_wrt_Wfom_1xN==trans(sum(dummy_mat_cx_NK_conj, 1));

			if (m_index==0){

				m_index_other=m_index_array[1];
			}

			if (m_index==1){

				m_index_other=m_index_array[0];
			}

			Pwrt_C_Wfom_fom(f_index, o_index, m_index)=real(accu( -G_ba_FxN_cx_mat.row(f_index)%(P_conj_Xhm_wrt_Wfom_1xN%conj((*Xhat_fnm_p).slice(m_index_other).row(f_index)))-conj(G_ba_FxN_cx_mat.row(f_index))%(P_Xhm_wrt_Wfom_1xN%((*Xhat_fnm_p).slice(m_index_other).row(f_index)))+2*(P_Xhm_wrt_Wfom_1xN%(*Xhat_fnm_p).slice(m_index).row(f_index)%abs((*Xhat_fnm_p).slice(m_index_other).row(f_index))) ));

			}	

		}	

}

/*} */

}

void complex_argument_costfun_m2_W_update3_entry(arg_struct_t* argStruct_p){

compute_Pwrt_C_Wfom_fom(1, 0, argStruct_p->expj_Phi_W_fom_p, argStruct_p->Y_lk_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p, argStruct_p->Z_ol_p, argStruct_p->Xhat_fnm_p, argStruct_p->Xtilde_fnm_p, argStruct_p->expj_Phi_S_nkf_p);

W_update3_Wfom(argStruct_p->W_fom_p);

}