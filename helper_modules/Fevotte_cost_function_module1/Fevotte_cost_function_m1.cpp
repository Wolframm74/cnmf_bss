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

static cx_cube::fixed<F_static, N_static, L_static> Y_est_FxNxL;
static cx_mat::fixed<F_static, N_static> entry_point1_fun2_mat1_FN;

static void entry_point1_fun2_compute_Y_est(mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p){

int l_index, k_index;

for (l_index=0; l_index<L_static; l_index++){

	entry_point1_fun2_mat1_FN.zeros();

	for (k_index=0; k_index<K_static; k_index++){

		entry_point1_fun2_mat1_FN=entry_point1_fun2_mat1_FN+((*Y_lk_p)(l_index, k_index))*expj_Phi_S_fnk_local.slice(k_index)%kron((*T_fk_p).col(k_index), trans((*V_nk_p).col(k_index)));

	}


Y_est_FxNxL.slice(l_index)=entry_point1_fun2_mat1_FN;

}

}

static cx_cube::fixed<M_static, F_static, N_static> Xtilde_mfn_local;
static cx_mat::fixed<F_static*N_static, M_static> entry_point1_fun3_fun0_mat1_FNxM;
static cx_mat::fixed<M_static, F_static*N_static> entry_point1_fun3_fun0_mat2_MxFN;

static cx_cube::fixed<L_static, F_static, N_static> Y_est_LxFxN_local;
static cx_mat::fixed<F_static*N_static, L_static> entry_point1_fun3_fun0_mat3_FNxL;
static cx_mat::fixed<L_static, F_static*N_static> entry_point1_fun3_fun0_mat4_LxFN;


static void entry_point1_fun3_fun0_rotate_X_Y(cx_cube* Xtilde_fnm_p){

/*copy the nkf input into NKxF */
arrayops::copy(entry_point1_fun3_fun0_mat1_FNxM.memptr(), (*Xtilde_fnm_p).memptr(), F_static*N_static*M_static);

/*transpose NKxF into FxNK*/	
entry_point1_fun3_fun0_mat2_MxFN=strans(entry_point1_fun3_fun0_mat1_FNxM);

/*copy FxNK into output FNK*/	
arrayops::copy(Xtilde_mfn_local.memptr(), entry_point1_fun3_fun0_mat2_MxFN.memptr(), F_static*N_static*M_static);

/*copy the nkf input into NKxF */
arrayops::copy(entry_point1_fun3_fun0_mat3_FNxL.memptr(), Y_est_FxNxL.memptr(), F_static*N_static*L_static);

/*transpose NKxF into FxNK*/	
entry_point1_fun3_fun0_mat4_LxFN=strans(entry_point1_fun3_fun0_mat3_FNxL);

/*copy FxNK into output FNK*/	
arrayops::copy(Y_est_LxFxN_local.memptr(), entry_point1_fun3_fun0_mat4_LxFN.memptr(), F_static*N_static*L_static);


}

static cx_cube::fixed<M_static, L_static, F_static_admissible> A_pseudoinv_MLF; 
static cx_mat::fixed<M_static, M_static> entry_point1_fun3_fun1_Sigma_B_MM_f;

static cx_mat::fixed<M_static, L_static> entry_point1_fun3_fun1_Rhat_XS_ML_f_N_accumat;
static cx_mat::fixed<L_static, L_static> entry_point1_fun3_fun1_Rhat_SS_LL_f_N_accumat;
static cx_mat::fixed<M_static, M_static> entry_point1_fun3_fun1_Rhat_XX_MM_f_N_accumat;

static cx_mat::fixed<M_static, L_static> entry_point1_fun3_fun1_Rhat_XS_ML_f;
static cx_mat::fixed<L_static, L_static> entry_point1_fun3_fun1_Rhat_SS_LL_f;
static cx_mat::fixed<M_static, M_static> entry_point1_fun3_fun1_Rhat_XX_MM_f;

static cx_mat::fixed<M_static, M_static> entry_point1_fun3_fun1_Sigma_X_MM_fn;
static cx_mat::fixed<L_static, L_static> entry_point1_fun3_fun1_Sigma_S_LL_fn;


static double C1_Fevotte_cost_function_output_variable;
static colvec::fixed<3>C1_Fevotte_cost_function_output_variable_colvec;

static void entry_point1_fun3_fun1_compute_A(void){

int f_index, n_index; 

int f_index_admissible;

double val_logdet;
double sign_logdet;

C1_Fevotte_cost_function_output_variable=0; 
C1_Fevotte_cost_function_output_variable_colvec.zeros();

for (f_index=0; f_index<F_static_admissible; f_index++){

	f_index_admissible=f_index+1; 

	entry_point1_fun3_fun1_Rhat_XX_MM_f_N_accumat.zeros();
	entry_point1_fun3_fun1_Rhat_XS_ML_f_N_accumat.zeros();
	entry_point1_fun3_fun1_Rhat_SS_LL_f_N_accumat.zeros();

	for (n_index=0; n_index<N_static; n_index++){

	entry_point1_fun3_fun1_Rhat_XX_MM_f_N_accumat=entry_point1_fun3_fun1_Rhat_XX_MM_f_N_accumat+Xtilde_mfn_local.slice(n_index).col(f_index_admissible)*trans(Xtilde_mfn_local.slice(n_index).col(f_index_admissible));

	entry_point1_fun3_fun1_Rhat_XS_ML_f_N_accumat=entry_point1_fun3_fun1_Rhat_XS_ML_f_N_accumat+Xtilde_mfn_local.slice(n_index).col(f_index_admissible)*trans(Y_est_LxFxN_local.slice(n_index).col(f_index_admissible));

	entry_point1_fun3_fun1_Rhat_SS_LL_f_N_accumat=entry_point1_fun3_fun1_Rhat_SS_LL_f_N_accumat+Y_est_LxFxN_local.slice(n_index).col(f_index_admissible)*trans(Y_est_LxFxN_local.slice(n_index).col(f_index_admissible));

/*	entry_point1_fun3_fun1_Sigma_B_MM_fn.diag()=*/

	}

	entry_point1_fun3_fun1_Rhat_XX_MM_f=(1/((double)N_static))*entry_point1_fun3_fun1_Rhat_XX_MM_f_N_accumat;

	entry_point1_fun3_fun1_Rhat_XS_ML_f=(1/((double)N_static))*entry_point1_fun3_fun1_Rhat_XS_ML_f_N_accumat;

	entry_point1_fun3_fun1_Rhat_SS_LL_f=(1/((double)N_static))*entry_point1_fun3_fun1_Rhat_SS_LL_f_N_accumat;

	A_pseudoinv_MLF.slice(f_index)=entry_point1_fun3_fun1_Rhat_XS_ML_f*pinv(entry_point1_fun3_fun1_Rhat_SS_LL_f);

	entry_point1_fun3_fun1_Sigma_B_MM_f=diagmat(entry_point1_fun3_fun1_Rhat_XX_MM_f- (A_pseudoinv_MLF.slice(f_index))*trans(entry_point1_fun3_fun1_Rhat_XS_ML_f) - entry_point1_fun3_fun1_Rhat_XS_ML_f*trans(A_pseudoinv_MLF.slice(f_index)) + A_pseudoinv_MLF.slice(f_index)*entry_point1_fun3_fun1_Rhat_SS_LL_f*trans(A_pseudoinv_MLF.slice(f_index)));

	for (n_index=0; n_index<N_static; n_index++){

	entry_point1_fun3_fun1_Sigma_S_LL_fn=diagmat(Y_est_LxFxN_local.slice(n_index).col(f_index));

	entry_point1_fun3_fun1_Sigma_X_MM_fn=A_pseudoinv_MLF.slice(f_index)*entry_point1_fun3_fun1_Sigma_S_LL_fn*trans(A_pseudoinv_MLF.slice(f_index))+entry_point1_fun3_fun1_Sigma_B_MM_f;

	C1_Fevotte_cost_function_output_variable_colvec(1)=C1_Fevotte_cost_function_output_variable_colvec(1)+real(trace(Xtilde_mfn_local.slice(n_index).col(f_index_admissible)*trans(Xtilde_mfn_local.slice(n_index).col(f_index_admissible))*pinv(entry_point1_fun3_fun1_Sigma_X_MM_fn)));

	C1_Fevotte_cost_function_output_variable_colvec(2)=C1_Fevotte_cost_function_output_variable_colvec(2)+log(real(det(entry_point1_fun3_fun1_Sigma_X_MM_fn)));

	C1_Fevotte_cost_function_output_variable=C1_Fevotte_cost_function_output_variable+real(trace(Xtilde_mfn_local.slice(n_index).col(f_index_admissible)*trans(Xtilde_mfn_local.slice(n_index).col(f_index_admissible))*pinv(entry_point1_fun3_fun1_Sigma_X_MM_fn)))+log(real(det(entry_point1_fun3_fun1_Sigma_X_MM_fn)));

	}

}

C1_Fevotte_cost_function_output_variable_colvec(0)=C1_Fevotte_cost_function_output_variable;

C1_Fevotte_cost_function_output_variable_colvec.print("C1_Fevotte_cost_function_output_variable_colvec");

}

static void entry_point1_fun3_fun1_compute_sigma_S_sigma_X(void){


}

static void entry_point1_fun3_compute_A_and_cost(cx_cube* Xtilde_fnm_p){

entry_point1_fun3_fun0_rotate_X_Y(Xtilde_fnm_p);

entry_point1_fun3_fun1_compute_A();

}

void Fevotte_cost_function_module1_entry_point1(cx_cube* Xtilde_fnm_p, mat* Z_ol_p, mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, cube* W_fom_p, cx_cube* expj_Phi_W_fom_p, cx_cube* expj_Phi_S_nkf_p){

entry_point1_fun1_rotate_Phi_S(expj_Phi_S_nkf_p);

entry_point1_fun2_compute_Y_est(Y_lk_p, T_fk_p, V_nk_p);

entry_point1_fun3_compute_A_and_cost(Xtilde_fnm_p);
	
}