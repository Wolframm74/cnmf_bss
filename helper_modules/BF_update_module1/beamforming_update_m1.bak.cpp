#include "local_inc.hpp"

static cx_cube::fixed<F_static, N_static, K_static> expj_Phi_S_fnk_local;
static cx_mat::fixed<N_static*K_static, F_static> entry_point1_fun1_mat1_NKxF;
static cx_mat::fixed<F_static, N_static*K_static> entry_point1_fun1_mat2_FxNK;

static void entry_point1_fun1_rotate_Phi_S(cx_cube* expj_Phi_S_nkf_p){

/*copy the nkf input into NKxF */
arrayops::copy(entry_point1_fun1_mat1_NKxF.memptr(), (*expj_Phi_S_nkf_p).memptr(), F_static*N_static*K_static);

/*transpose NKxF into FxNK*/	
entry_point1_fun1_mat2_FxNK=strans(entry_point1_fun1_mat1_NKxF);

/*copy FxNK into output FNK*/	
arrayops::copy(expj_Phi_S_fnk_local.memptr(), entry_point1_fun1_mat2_FxNK.memptr(), F_static*N_static*K_static);

}

static rowvec::fixed<K_static> frob_norm_rowvec_norm_V;

static void normalize_V_nk_local(mat* V_nk_p){

frob_norm_rowvec_norm_V=sqrt(sum( (*V_nk_p)%(*V_nk_p) , 0));

/*Divide out each col by its frob norm:*/
(*V_nk_p)=(*V_nk_p)/kron(ones_col_Nx1, frob_norm_rowvec_norm_V);

}

static void mat::fixed<N_static, K_static> V_nk_complement; 

static void entry_point1_fun2_compute_V_complement(mat* V_nk_p, mat* T_fk_p){

mat* V_nk_complement_p=&V_nk_complement;

/*first normalize V*/
normalize_V_nk_frob_norm(V_nk_p, T_fk_p);

/*Set the complement to ones*/

V_nk_complement.ones();

/*Subtract off V from the complement*/
V_nk_complement=V_nk_complement-(*V_nk_p);

/*There should only be ones left in positions where V was originally zero. This is what you want*/

/*Normalize the complement to bring the average energy down a bit*/
normalize_V_nk_local(V_nk_complement_p);

}


static fixed::colvec<O_static> kith_dummy_sum_vec;
static fixed::mat<K_static, O_static> Z_ko_target_BF;

entry_point1_fun3_compute_Zko_target(mat* T_fk_p, mat* V_nk_p, cube* W_fom_p, cx_cube* expj_Phi_W_fom_p, cx_cube* expj_Phi_S_nkf_p){

Z_ko_target_BF.zeros();

int k_i_index, k_index;

int m_index, f_index, n_index, o_index; 

uword min_o_index;

for (k_i_index=0; k_i_index<K_static; k_i_index++){

	kith_dummy_sum_vec.zeros();

	for (o_index=0; o_index<O_static; o_index++){



		for (m_index=0; m_index<M_static; m_index++){

			for (f_index=0; f_index<F_static; f_index++){

				for (n_index=0; n_index<N_static; n_index++){

					for (k_index=0; k_index<K_static; k_index++){

						if (k_index!=k_i_index){

						kith_dummy_sum_vec(o_index)=kith_dummy_sum_vec(o_index)+pow( abs ( ((*W_fom_p)(f_index, o_index, m_index))*((*expj_Phi_W_fom_p)(f_index, o_index, m_index))*((*T_fk_p)(f_index, k_index))*(V_nk_complement(n_index, k_index))*((*V_nk_p)(n_index, k_index))((*expj_Phi_S_nkf_p)(n_index, k_index, f_index)) ) , 2);

						}

					}

				}

			}

		}

	}

	kith_dummy_sum_vec.min(min_o_index);

	Z_ko_target_BF(k_i_index, (int)min_o_index)=1;

}


}

void Beamforming_update_module1_entry_point1(mat* T_fk_p, mat* V_nk_p, cube* W_fom_p, cx_cube* expj_Phi_W_fom_p, cx_cube* expj_Phi_S_nkf_p){

/*Rotate Phi_S*/
/*entry_point1_fun1_rotate_Phi_S(expj_Phi_S_nkf_p);*/

/*Rotate function might be useful if you wanted to create a code that used accu(FxN) instead of so many for loops as done above. */	

/*Compute the complement of V*/	
entry_point1_fun2_compute_V_complement(V_nk_p, T_fk_p);

/*Iterate for k=1:K, o=1:O*/	
entry_point1_fun3_compute_Zko_target(mat* T_fk_p, mat* V_nk_p, cube* W_fom_p, cx_cube* expj_Phi_W_fom_p, cx_cube* expj_Phi_S_nkf_p);

/*This should compute a target matrix Z_ko target*/	

/*Integrate Z_ko target into the current model*/	

}