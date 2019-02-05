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



/*void Fevotte_cost_function_module1_entry_point1(cx_cube* Xtilde_fnm_p, mat* Z_ol_p, mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, cube* W_fom_p, cx_cube* expj_Phi_W_fom_p, cx_cube* expj_Phi_S_nkf_p){

entry_point1_fun1_rotate_Phi_S(expj_Phi_S_nkf_p);

entry_point1_fun2_compute_Y_est(Y_lk_p, T_fk_p, V_nk_p);

entry_point1_fun3_compute_A_and_cost(Xtilde_fnm_p);
	
}*/

static cube::fixed<F_static, N_static, K_static> U_fnk;

static void module2_compute_U_fnk(mat* T_fk_p, mat* V_nk_p){

int k_index;

for (k_index=0; k_index<K_static; k_index++){

	U_fnk.slice(k_index)=square(kron((*T_fk_p).col(k_index), trans((*V_nk_p).col(k_index))));

}

}

static void module2_T_fk_update(mat* T_fk_p, mat* V_nk_p){

int k_index, f_index;

module2_compute_U_fnk(T_fk_p, V_nk_p);

for (f_index=0; f_index<F_static; f_index++){

	for (k_index=0; k_index<K_static; k_index++){

		/*(*T_fk_p)(f_index, k_index)=0.5*(*T_fk_p)(f_index, k_index)+0.5*(1/((double)N_static))*as_scalar(sum( trans(U_fnk.slice(k_index).row(f_index)) /( (*V_nk_p).col(k_index) ) ));*/
		(*T_fk_p)(f_index, k_index)=0.8*(*T_fk_p)(f_index, k_index)+0.2*(1/((double)N_static))*as_scalar(sum( trans(U_fnk.slice(k_index).row(f_index)) /( (*V_nk_p).col(k_index) ) ));

	}

}

}

static void module2_V_nk_update(mat* T_fk_p, mat* V_nk_p){

int k_index, n_index;

module2_compute_U_fnk(T_fk_p, V_nk_p);

for (n_index=0; n_index<N_static; n_index++){

	for (k_index=0; k_index<K_static; k_index++){

		(*V_nk_p)(n_index, k_index)=0.8*(*V_nk_p)(n_index, k_index)+0.2*(1/((double)F_static))*as_scalar(sum( (U_fnk.slice(k_index).col(n_index)) /( (*T_fk_p).col(k_index) ) ));

	}

}

}

void Fevotte_cost_function_module2_T_fk_update(mat* T_fk_p, mat* V_nk_p, cx_cube* expj_Phi_S_nkf_p){

entry_point1_fun1_rotate_Phi_S(expj_Phi_S_nkf_p);

module2_T_fk_update(T_fk_p, V_nk_p);
	
}

void Fevotte_cost_function_module2_V_nk_update(mat* T_fk_p, mat* V_nk_p, cx_cube* expj_Phi_S_nkf_p){

entry_point1_fun1_rotate_Phi_S(expj_Phi_S_nkf_p);

module2_V_nk_update(T_fk_p, V_nk_p);
	
}