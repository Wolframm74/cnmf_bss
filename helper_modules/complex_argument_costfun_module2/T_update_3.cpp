#include "local_inc.hpp"

cx_cube::fixed<F_static, K_static, N_static> P_Xh1_wrt_Tfk_fkn;

static void populate_P_Xh1_wrt_Tfk_fkn(int m_index, cx_cube* W_fom_cx_p, mat* Z_ol_p, mat* Y_lk_p, cx_cube* expj_Phi_S_fkn_p, mat* T_fk_p, mat* V_nk_p){

int k_index, n_index; 

populate_dummy_mat_H_fk_1(m_index, W_fom_cx_p, Z_ol_p, Y_lk_p);

for (n_index=0; n_index<N_static; n_index++){

	for (k_index=0; k_index<K_static; k_index++){

	P_Xh1_wrt_Tfk_fkn.slice(n_index).col(k_index)=((*V_nk_p)(n_index, k_index))*(dummy_mat_H_fk_1.col(k_index))%((*expj_Phi_S_fkn_p).slice(n_index).col(k_index));

	}

}

}


/*Function local data*/

/*

The following functions and their data used to reside here, but now moved to "common_module.cpp"

static void populate_G_H_fn_channel1(int m_index, cx_cube* Xhat_fnm_p);

static void populate_G_H_fn_channel2(int m_index, cx_cube* Xhat_fnm_p);

*/

/*Function local data*/
mat::fixed<F_static, N_static> T_update_dummy_mat_real_FN;

/*Function local data*/

cube::fixed<F_static, K_static, N_static> P_G1_wrt_Tfk_fkn;

static void populate_P_G1_wrt_Tfk_fkn(void){

/*Take the imaginary part of Pwrt_Xh1_fnk*/
P_G1_wrt_Tfk_fkn=imag(P_Xh1_wrt_Tfk_fkn);


}

/*Function local data*/

cube::fixed<F_static, K_static, N_static> P_H1_wrt_Tfk_fkn;

cube::fixed<F_static, K_static, N_static> P_H1_wrt_Tfk_num_fkn;

static void populate_P_H1_wrt_Tfk_fkn(int m_index, cx_cube* Xhat_fnm_p){

int n_index, k_index;

T_update_dummy_mat_real_FN=real((*Xhat_fnm_p).slice(m_index));

for (n_index=0; n_index<N_static; n_index++){

	for (k_index=0; k_index<K_static; k_index++){

	P_H1_wrt_Tfk_num_fkn.slice(n_index).col(k_index)=(T_update_dummy_mat_real_FN.col(n_index))%(real(P_Xh1_wrt_Tfk_fkn.slice(n_index).col(k_index)));

	P_H1_wrt_Tfk_num_fkn.slice(n_index).col(k_index)=(G1_outmat_fn.col(n_index))%(imag(P_Xh1_wrt_Tfk_fkn.slice(n_index).col(k_index)))+(P_H1_wrt_Tfk_num_fkn.slice(n_index).col(k_index));

	P_H1_wrt_Tfk_fkn.slice(n_index).col(k_index)=(P_H1_wrt_Tfk_num_fkn.slice(n_index).col(k_index))/(sqrt_denmat_1_fn.col(n_index));

	}
	
}

}

/*Function local data*/

mat::fixed<F_static, K_static> dummy_mat1_real_FK;	/*Store the G1, G2 related stuff*/
mat::fixed<F_static, K_static> dummy_mat2_real_FK;	/*Store the H1, H2 related stuff*/

cube::fixed<F_static, K_static, N_static> P_arg_Xh1_wrt_Tfk_fkn;

static void populate_P_arg_Xh1_wrt_Tfk_fkn(void){

int n_index;

for (n_index=0; n_index<N_static; n_index++){

dummy_mat1_real_FK=kron(G1_outmat_fn.col(n_index), ones_row_1xK);

dummy_mat2_real_FK=kron(H1_outmat_fn.col(n_index), ones_row_1xK);

P_arg_Xh1_wrt_Tfk_fkn.slice(n_index)=2*ones_cube_FK/(square(dummy_mat1_real_FK)+square(dummy_mat2_real_FK));

P_arg_Xh1_wrt_Tfk_fkn.slice(n_index)=(P_arg_Xh1_wrt_Tfk_fkn.slice(n_index))%(((P_G1_wrt_Tfk_fkn.slice(n_index))%dummy_mat2_real_FK)-(dummy_mat1_real_FK%(P_H1_wrt_Tfk_fkn.slice(n_index))));

}

}

/*Function local data*/
cube::fixed<F_static, K_static, N_static> Pwrt_C_Tfk_fkn;

mat::fixed<F_static, K_static> Pwrt_C_Tfk_fk_mth_component;

mat::fixed<F_static, K_static> Pwrt_C_Tfk_fk_outmat;

static mat& compute_Pwrt_C_Tfk_fkn(int m_index, cx_cube* Xhat_fnm_p, cx_cube* Xtilde_fnm_p){

int n_index;

mat& Pwrt_C_Tfk_fk_mth_component_ref=Pwrt_C_Tfk_fk_mth_component;

Pwrt_C_Tfk_fk_mth_component.zeros();

/*Elementise wise FxNxK product of the thing you computed with the previous function vs the difference between Pwrt_arg_Xh1_fnk and Pwrt_arg_Xh2_fnk */
T_update_dummy_mat_real_FN.zeros();

T_update_dummy_mat_real_FN=compute_arg_X_FN(Xhat_fnm_p, m_index);

/*T_update_dummy_mat_real_FN=T_update_dummy_mat_real_FN+compute_arg_X_FN(Xtilde_fnm_p, m_index_2);*/

T_update_dummy_mat_real_FN=T_update_dummy_mat_real_FN-compute_arg_X_FN(Xtilde_fnm_p, m_index);

/*T_update_dummy_mat_real_FN=T_update_dummy_mat_real_FN-compute_arg_X_FN(Xhat_fnm_p, m_index_2);*/

T_update_dummy_mat_real_FN=2*T_update_dummy_mat_real_FN;

/*Spead the result over k=0:K-1*/
for (n_index=0; n_index<N_static; n_index++){

dummy_mat1_real_FK=kron(T_update_dummy_mat_real_FN.col(n_index), ones_row_1xK);

Pwrt_C_Tfk_fkn.slice(n_index)=dummy_mat1_real_FK%(P_arg_Xh1_wrt_Tfk_fkn.slice(n_index));

Pwrt_C_Tfk_fk_mth_component=Pwrt_C_Tfk_fk_mth_component+Pwrt_C_Tfk_fkn.slice(n_index);

}	

return Pwrt_C_Tfk_fk_mth_component_ref;

}

static void T_update2_Tfk(mat* T_fk_p){

double epsilon_Tfk;

epsilon_Tfk=0.01;

(*T_fk_p)=(*T_fk_p)-(epsilon_Tfk)*(Pwrt_C_Tfk_fk_outmat);

}

void T_update2_entry(arg_struct_t* argStruct_p){

int m_index;
/*int m_index_2=1;*/

Pwrt_C_Tfk_fk_outmat.zeros();

for (m_index=0; m_index<M_static; m_index++){

populate_P_Xh1_wrt_Tfk_fkn(m_index, argStruct_p->W_fom_cx_p, argStruct_p->Z_ol_p, argStruct_p->Y_lk_p, argStruct_p->expj_Phi_S_fkn_p, argStruct_p->T_fk_p, argStruct_p->V_nk_p);

populate_G_H_fn_channel1(m_index, argStruct_p->Xhat_fnm_p);

populate_P_G1_wrt_Tfk_fkn();

populate_P_H1_wrt_Tfk_fkn(m_index, argStruct_p->Xhat_fnm_p);

populate_P_arg_Xh1_wrt_Tfk_fkn();

Pwrt_C_Tfk_fk_outmat=Pwrt_C_Tfk_fk_outmat+compute_Pwrt_C_Tfk_fkn(m_index, argStruct_p->Xhat_fnm_p, argStruct_p->Xtilde_fnm_p);

}

T_update2_Tfk(argStruct_p->T_fk_p);

}

















/*Would like the output of everything to perhaps be in R: FxKxN because you eventually need to integrate out n=0:N-1 */

/*To spreak the arg(X_fn)'s across k should be easy, simply index over n=0:N-1*/

/*To populate the P_Xh1_wrt_T tensor , simply do a double for loop over k, n, to index out the Fx1 column vectors*/
