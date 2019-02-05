#include "local_inc.hpp"

#define NUM_FEVOTTE_EM_ITERATIONS 2

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

static cx_cube::fixed<F_static, N_static, K_static> Y_est_FxNxK;
static cx_cube::fixed<F_static, N_static, L_static> Y_est_FxNxL;
static cx_mat::fixed<F_static, N_static> entry_point1_fun2_mat1_FN;

static cube::fixed<F_static, N_static, K_static> Y_est_sigma_input_FxNxK;
static cube::fixed<F_static, N_static, L_static> Y_est_sigma_input_FxNxL;
static mat::fixed<F_static, N_static> entry_point1_fun2_mat2_FN;

static void entry_point1_fun2_compute_Y_est(mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p){

int l_index, k_index;

for (l_index=0; l_index<L_static; l_index++){

	entry_point1_fun2_mat1_FN.zeros();

	for (k_index=0; k_index<K_static; k_index++){

		entry_point1_fun2_mat1_FN=entry_point1_fun2_mat1_FN+((*Y_lk_p)(l_index, k_index))*expj_Phi_S_fnk_local.slice(k_index)%kron((*T_fk_p).col(k_index), trans((*V_nk_p).col(k_index)));
		entry_point1_fun2_mat2_FN=entry_point1_fun2_mat2_FN+((*Y_lk_p)(l_index, k_index))*kron((*T_fk_p).col(k_index), trans((*V_nk_p).col(k_index)));

	}


Y_est_FxNxL.slice(l_index)=entry_point1_fun2_mat1_FN;
Y_est_sigma_input_FxNxL.slice(l_index)=entry_point1_fun2_mat2_FN;

}

for (k_index=0; k_index<K_static; k_index++){

	Y_est_FxNxK.slice(k_index)=expj_Phi_S_fnk_local.slice(k_index)%kron((*T_fk_p).col(k_index), trans((*V_nk_p).col(k_index)));
	Y_est_sigma_input_FxNxK.slice(k_index)=kron((*T_fk_p).col(k_index), trans((*V_nk_p).col(k_index)));

}

/*Y_est_FxNxL.print("line 58");

Y_est_FxNxK.print("line 60");*/

}

static cx_cube::fixed<M_static, F_static, N_static> Xtilde_mfn_local;
static cx_mat::fixed<F_static*N_static, M_static> entry_point1_fun0_mat1_FNxM;
static cx_mat::fixed<M_static, F_static*N_static> entry_point1_fun0_mat2_MxFN;

static cx_cube::fixed<L_static, F_static, N_static> Y_est_LxFxN_local;
static cx_mat::fixed<F_static*N_static, L_static> entry_point1_fun0_mat3_FNxL;
static cx_mat::fixed<L_static, F_static*N_static> entry_point1_fun0_mat4_LxFN;

static cx_cube::fixed<K_static, F_static, N_static> Y_est_KxFxN_local;
static cx_mat::fixed<F_static*N_static, K_static> entry_point1_fun0_mat5_FNxK;
static cx_mat::fixed<K_static, F_static*N_static> entry_point1_fun0_mat6_KxFN;

static cube::fixed<L_static, F_static, N_static> Y_est_LxFxN_sigma_input_local;
static mat::fixed<F_static*N_static, L_static> entry_point1_fun0_mat7_FNxL;
static mat::fixed<L_static, F_static*N_static> entry_point1_fun0_mat8_LxFN;

static cube::fixed<K_static, F_static, N_static> Y_est_KxFxN_sigma_input_local;
static mat::fixed<F_static*N_static, K_static> entry_point1_fun0_mat9_FNxK;
static mat::fixed<K_static, F_static*N_static> entry_point1_fun0_mat10_KxFN;

static void entry_point1_fun0_rotate_X_Y(cx_cube* Xtilde_fnm_p){

/*copy the nkf input into NKxF */
arrayops::copy(entry_point1_fun0_mat1_FNxM.memptr(), (*Xtilde_fnm_p).memptr(), F_static*N_static*M_static);

/*transpose NKxF into FxNK*/	
entry_point1_fun0_mat2_MxFN=strans(entry_point1_fun0_mat1_FNxM);

/*copy FxNK into output FNK*/	
arrayops::copy(Xtilde_mfn_local.memptr(), entry_point1_fun0_mat2_MxFN.memptr(), F_static*N_static*M_static);

/*Y_est_FxNxL.slice(0).print("Y_est_FxNxL.slice(0), line 95; RIGHT BEFORE THE FIRST MEMCPY");*/

/*copy the nkf input into NKxF */
arrayops::copy(entry_point1_fun0_mat3_FNxL.memptr(), Y_est_FxNxL.memptr(), F_static*N_static*L_static);

/*transpose NKxF into FxNK*/	
entry_point1_fun0_mat4_LxFN=strans(entry_point1_fun0_mat3_FNxL);

/*copy FxNK into output FNK*/	
arrayops::copy(Y_est_LxFxN_local.memptr(), entry_point1_fun0_mat4_LxFN.memptr(), F_static*N_static*L_static);

/*Y_est_LxFxN_local.slice(0).print("Y_est_LxFxN_local.slice(0), line 104; RIGHT AFTER THE LAST MEMCPY");*/

/*copy the nkf input into NKxF */
arrayops::copy(entry_point1_fun0_mat5_FNxK.memptr(), Y_est_FxNxK.memptr(), F_static*N_static*K_static);

/*transpose NKxF into FxNK*/	
entry_point1_fun0_mat6_KxFN=strans(entry_point1_fun0_mat5_FNxK);

/*copy FxNK into output FNK*/	
arrayops::copy(Y_est_KxFxN_local.memptr(), entry_point1_fun0_mat6_KxFN.memptr(), F_static*N_static*K_static);


arrayops::copy(entry_point1_fun0_mat7_FNxL.memptr(), Y_est_sigma_input_FxNxL.memptr(), F_static*N_static*L_static);
entry_point1_fun0_mat8_LxFN=strans(entry_point1_fun0_mat7_FNxL);
arrayops::copy(Y_est_LxFxN_sigma_input_local.memptr(), entry_point1_fun0_mat8_LxFN.memptr(), F_static*N_static*L_static);

arrayops::copy(entry_point1_fun0_mat9_FNxK.memptr(), Y_est_sigma_input_FxNxK.memptr(), F_static*N_static*K_static);
entry_point1_fun0_mat10_KxFN=strans(entry_point1_fun0_mat9_FNxK);
arrayops::copy(Y_est_KxFxN_sigma_input_local.memptr(), entry_point1_fun0_mat10_KxFN.memptr(), F_static*N_static*K_static);

}

static cx_cube::fixed<M_static, L_static, F_static_admissible> A_pseudoinv_EM_MLF; 
static cx_cube::fixed<M_static, K_static, F_static_admissible> A_pseudoinv_EM_augmented_MKF; 
static cx_mat::fixed<M_static, M_static> entry_point1_fun3_Sigma_B_MM_f;

static cx_mat::fixed<M_static, L_static> entry_point1_fun3_Rhat_XS_ML_f_N_accumat;
static cx_mat::fixed<L_static, L_static> entry_point1_fun3_Rhat_SS_LL_f_N_accumat;
static cx_mat::fixed<M_static, M_static> entry_point1_fun3_Rhat_XX_MM_f_N_accumat;

static cx_mat::fixed<M_static, L_static> entry_point1_fun3_Rhat_XS_ML_f;
static cx_mat::fixed<L_static, L_static> entry_point1_fun3_Rhat_SS_LL_f;
static cx_mat::fixed<M_static, M_static> entry_point1_fun3_Rhat_XX_MM_f;

static cx_mat::fixed<M_static, M_static> entry_point1_fun3_Sigma_X_MM_fn;
static mat::fixed<L_static, L_static> entry_point1_fun3_Sigma_S_LL_fn;
static mat::fixed<K_static, K_static> entry_point1_fun3_Sigma_C_KK_fn;

static cx_mat::fixed<K_static, M_static> entry_point1_fun3_G_wiener_C_KM_fn;
static cx_mat::fixed<L_static, M_static> entry_point1_fun3_G_wiener_S_LM_fn;

static cx_cube::fixed<K_static, F_static_admissible, N_static> entry_point1_fun3_Chat_colvecs_K_fn;
static cx_cube::fixed<L_static, F_static_admissible, N_static> entry_point1_fun3_Shat_colvecs_L_fn;

static cx_cube::fixed<K_static, F_static_admissible, N_static> entry_point1_fun3_output_phase_Phi_S_K_fn;

static cube::fixed<K_static, F_static_admissible, N_static> entry_point1_fun3_uhat_KFN;
static mat::fixed<K_static, K_static> entry_point1_fun3_uhat_KK_fn;

static double C1_Fevotte_cost_function_output_variable;
static colvec::fixed<3>C1_Fevotte_cost_function_output_variable_colvec;

static colvec::fixed<L_static> entry_point1_fun3_kth_column_Y_lk; 

static void entry_point1_fun3_compute_A_via_EM(mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, cx_cube* expj_Phi_S_nkf_p, cx_cube* expj_Phi_S_fkn_p){

int f_index, n_index, k_index; 

uword l_index_mapping;

int f_index_admissible;

double val_logdet;
double sign_logdet;

bool init_flag=true; /*Flag set to true signifying A needs to be initialized*/

int em_counter; 

double dummy_accu_value_fk, dummy_accu_value_kn; 

C1_Fevotte_cost_function_output_variable=0; 
C1_Fevotte_cost_function_output_variable_colvec.zeros();

/*(*Y_lk_p).print("Y_lk_p, line 172");*/

/*Y_est_LxFxN_local.slice(0).print("Y_est_LxFxN_local.slice(0), line 178;");*/

for (em_counter=0; em_counter<NUM_FEVOTTE_EM_ITERATIONS; em_counter++){

	if (init_flag==false){

			/*Signma C and S must be recomputed given the newly computed T_fk and V_nk parameters. */

			entry_point1_fun2_compute_Y_est(Y_lk_p, T_fk_p, V_nk_p);

			/*Compute the rotated quantities */
			arrayops::copy(entry_point1_fun0_mat7_FNxL.memptr(), Y_est_sigma_input_FxNxL.memptr(), F_static*N_static*L_static);
			entry_point1_fun0_mat8_LxFN=strans(entry_point1_fun0_mat7_FNxL);
			arrayops::copy(Y_est_LxFxN_sigma_input_local.memptr(), entry_point1_fun0_mat8_LxFN.memptr(), F_static*N_static*L_static);

			arrayops::copy(entry_point1_fun0_mat9_FNxK.memptr(), Y_est_sigma_input_FxNxK.memptr(), F_static*N_static*K_static);
			entry_point1_fun0_mat10_KxFN=strans(entry_point1_fun0_mat9_FNxK);
			arrayops::copy(Y_est_KxFxN_sigma_input_local.memptr(), entry_point1_fun0_mat10_KxFN.memptr(), F_static*N_static*K_static);

	}

	for (f_index=0; f_index<F_static_admissible; f_index++){

		f_index_admissible=f_index+1;

				/*C1_Fevotte_cost_function_output_variable_colvec.print("Fevotte EM line 175: ARGHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111");*/

		if (em_counter==0){

			/*Estimate A_f crudely without yet knowledge of the exact statistics*/
			entry_point1_fun3_Rhat_XS_ML_f_N_accumat.zeros();
			entry_point1_fun3_Rhat_SS_LL_f_N_accumat.zeros();

			for (n_index=0; n_index<N_static; n_index++){

			/*Xtilde_mfn_local.slice(n_index).col(f_index_admissible).print("line 207:");
			Y_est_LxFxN_local.slice(n_index).col(f_index_admissible).print("line 208:");*/

			entry_point1_fun3_Rhat_XS_ML_f_N_accumat=entry_point1_fun3_Rhat_XS_ML_f_N_accumat+Xtilde_mfn_local.slice(n_index).col(f_index_admissible)*trans(Y_est_LxFxN_local.slice(n_index).col(f_index_admissible));

			entry_point1_fun3_Rhat_SS_LL_f_N_accumat=entry_point1_fun3_Rhat_SS_LL_f_N_accumat+Y_est_LxFxN_local.slice(n_index).col(f_index_admissible)*trans(Y_est_LxFxN_local.slice(n_index).col(f_index_admissible));

		/*	entry_point1_fun3_Sigma_B_MM_fn.diag()=*/

			}

			/*entry_point1_fun3_Rhat_XS_ML_f_N_accumat.print("entry_point1_fun3_Rhat_XS_ML_f_N_accumat, line 215:");
			entry_point1_fun3_Rhat_XS_ML_f_N_accumat.print("entry_point1_fun3_Rhat_XS_ML_f_N_accumat, line 216:");*/

			entry_point1_fun3_Rhat_XS_ML_f=(1/((double)N_static))*entry_point1_fun3_Rhat_XS_ML_f_N_accumat;

			entry_point1_fun3_Rhat_SS_LL_f=(1/((double)N_static))*entry_point1_fun3_Rhat_SS_LL_f_N_accumat;

			A_pseudoinv_EM_MLF.slice(f_index)=entry_point1_fun3_Rhat_XS_ML_f*pinv(entry_point1_fun3_Rhat_SS_LL_f);

			init_flag=false;

		}

		/*C1_Fevotte_cost_function_output_variable_colvec.print("Fevotte EM line 202: ARGHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111");*/

		for (k_index=0; k_index<K_static; k_index++){

		entry_point1_fun3_kth_column_Y_lk=(*Y_lk_p).col(k_index);

		entry_point1_fun3_kth_column_Y_lk.max(l_index_mapping);
		 
		A_pseudoinv_EM_augmented_MKF.slice(f_index).col(k_index)=A_pseudoinv_EM_MLF.slice(f_index).col((int)l_index_mapping); 

		}

		/*Compute the E-step. Assume that you now have a crude estimate of A_pseudoinv_EM_MLF, that should improve once the statistics are all computed*/

		f_index_admissible=f_index+1; 

		entry_point1_fun3_Rhat_XX_MM_f_N_accumat.zeros();
		entry_point1_fun3_Rhat_XS_ML_f_N_accumat.zeros();
		entry_point1_fun3_Rhat_SS_LL_f_N_accumat.zeros();

		/*A_pseudoinv_EM_MLF.slice(f_index).print("A_pseudoinv_EM_MLF.slice(f_index), line 245:");
		A_pseudoinv_EM_augmented_MKF.slice(f_index).print("A_pseudoinv_EM_augmented_MKF.slice(f_index), line 246:");*/

		for (n_index=0; n_index<N_static; n_index++){

		if (init_flag==true){

			entry_point1_fun3_Sigma_X_MM_fn=A_pseudoinv_EM_MLF.slice(f_index)*entry_point1_fun3_Sigma_S_LL_fn*trans(A_pseudoinv_EM_MLF.slice(f_index));	/*Optionnally add the time independent, frequency dependent b noise term if you're returning from an M step*/

			entry_point1_fun3_Sigma_C_KK_fn=diagmat(Y_est_KxFxN_sigma_input_local.slice(n_index).col(f_index));

			entry_point1_fun3_Sigma_S_LL_fn=diagmat(Y_est_LxFxN_sigma_input_local.slice(n_index).col(f_index));

			/*Covariance matrices must be Hermitian symmetric, ie: real valued along main diagonal...*/

		} else if (init_flag==false){

			entry_point1_fun3_Sigma_X_MM_fn=A_pseudoinv_EM_MLF.slice(f_index)*entry_point1_fun3_Sigma_S_LL_fn*trans(A_pseudoinv_EM_MLF.slice(f_index))+entry_point1_fun3_Sigma_B_MM_f;

			/*Can now compute the outputs similarly as in the other case. */
			entry_point1_fun3_Sigma_C_KK_fn=diagmat(Y_est_KxFxN_sigma_input_local.slice(n_index).col(f_index));

			entry_point1_fun3_Sigma_S_LL_fn=diagmat(Y_est_LxFxN_sigma_input_local.slice(n_index).col(f_index));

			/*Covariance matrices must be Hermitian symmetric, ie: real valued along main diagonal...*/

		}			

		/*C1_Fevotte_cost_function_output_variable_colvec.print("Fevotte EM line 258: ARGHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111");*/

		entry_point1_fun3_G_wiener_S_LM_fn=entry_point1_fun3_Sigma_S_LL_fn*trans(A_pseudoinv_EM_MLF.slice(f_index))*pinv(entry_point1_fun3_Sigma_X_MM_fn);

		entry_point1_fun3_G_wiener_C_KM_fn=entry_point1_fun3_Sigma_C_KK_fn*trans(A_pseudoinv_EM_augmented_MKF.slice(f_index))*pinv(entry_point1_fun3_Sigma_X_MM_fn);
		
		entry_point1_fun3_Chat_colvecs_K_fn.slice(n_index).col(f_index)=entry_point1_fun3_G_wiener_C_KM_fn*Xtilde_mfn_local.slice(n_index).col(f_index_admissible);
		entry_point1_fun3_Shat_colvecs_L_fn.slice(n_index).col(f_index)=entry_point1_fun3_G_wiener_S_LM_fn*Xtilde_mfn_local.slice(n_index).col(f_index_admissible);

		/*entry_point1_fun3_Chat_colvecs_K_fn.slice(n_index).col(f_index).print("line 277: entry_point1_fun3_Chat_colvecs_K_fn.slice(n_index).col(f_index)");
		entry_point1_fun3_Shat_colvecs_L_fn.slice(n_index).col(f_index).print("line 278: entry_point1_fun3_Shat_colvecs_L_fn.slice(n_index).col(f_index)");*/

		entry_point1_fun3_Rhat_XX_MM_f_N_accumat=entry_point1_fun3_Rhat_XX_MM_f_N_accumat+Xtilde_mfn_local.slice(n_index).col(f_index_admissible)*trans(Xtilde_mfn_local.slice(n_index).col(f_index_admissible));

		entry_point1_fun3_Rhat_XS_ML_f_N_accumat=entry_point1_fun3_Rhat_XS_ML_f_N_accumat+Xtilde_mfn_local.slice(n_index).col(f_index_admissible)*trans(entry_point1_fun3_Shat_colvecs_L_fn.slice(n_index).col(f_index));

		entry_point1_fun3_Rhat_SS_LL_f_N_accumat=entry_point1_fun3_Rhat_SS_LL_f_N_accumat+entry_point1_fun3_Shat_colvecs_L_fn.slice(n_index).col(f_index)*trans(entry_point1_fun3_Shat_colvecs_L_fn.slice(n_index).col(f_index))+entry_point1_fun3_Sigma_S_LL_fn-(entry_point1_fun3_G_wiener_S_LM_fn*A_pseudoinv_EM_MLF.slice(f_index)*entry_point1_fun3_Sigma_S_LL_fn);

		entry_point1_fun3_uhat_KK_fn=real(entry_point1_fun3_Chat_colvecs_K_fn.slice(n_index).col(f_index)*trans(entry_point1_fun3_Chat_colvecs_K_fn.slice(n_index).col(f_index))+entry_point1_fun3_Sigma_C_KK_fn-(entry_point1_fun3_G_wiener_C_KM_fn*A_pseudoinv_EM_augmented_MKF.slice(f_index)*entry_point1_fun3_Sigma_C_KK_fn));

		/*entry_point1_fun3_uhat_KK_fn.print("entry_point1_fun3_uhat_KK_fn, line 285:");*/

		entry_point1_fun3_uhat_KFN.slice(n_index).col(f_index)=diagvec(entry_point1_fun3_uhat_KK_fn);

		/*entry_point1_fun3_uhat_KFN.print("entry_point1_fun3_uhat_KFN, line 289:");*/

	/*	entry_point1_fun3_Sigma_B_MM_fn.diag()=*/

		}	/*end-for: n_index*/

		entry_point1_fun3_Rhat_XX_MM_f=(1/((double)N_static))*entry_point1_fun3_Rhat_XX_MM_f_N_accumat;

		entry_point1_fun3_Rhat_XS_ML_f=(1/((double)N_static))*entry_point1_fun3_Rhat_XS_ML_f_N_accumat;

		entry_point1_fun3_Rhat_SS_LL_f=(1/((double)N_static))*entry_point1_fun3_Rhat_SS_LL_f_N_accumat;

		/*Compute the M-Step*/

		A_pseudoinv_EM_MLF.slice(f_index)=entry_point1_fun3_Rhat_XS_ML_f*pinv(entry_point1_fun3_Rhat_SS_LL_f);

		entry_point1_fun3_Sigma_B_MM_f=diagmat(entry_point1_fun3_Rhat_XX_MM_f- (A_pseudoinv_EM_MLF.slice(f_index))*trans(entry_point1_fun3_Rhat_XS_ML_f) - entry_point1_fun3_Rhat_XS_ML_f*trans(A_pseudoinv_EM_MLF.slice(f_index)) + A_pseudoinv_EM_MLF.slice(f_index)*entry_point1_fun3_Rhat_SS_LL_f*trans(A_pseudoinv_EM_MLF.slice(f_index)));

	}	/*end-for: f_index*/


/*(*T_fk_p).print("T_fk, line 306: before shot");*/


	/*update T_fk and V_kn as a function of k( for k=1:K) */
	for (k_index=0; k_index<K_static; k_index++){

		for (f_index=0; f_index<F_static_admissible; f_index++){

			f_index_admissible=f_index+1; 

			dummy_accu_value_fk=0;

			for (n_index=0; n_index<N_static; n_index++){

			dummy_accu_value_fk=dummy_accu_value_fk+entry_point1_fun3_uhat_KFN(k_index, f_index, n_index)/((*V_nk_p)(n_index, k_index));

			}

			dummy_accu_value_fk=(1/((double)N_static))*dummy_accu_value_fk;

			(*T_fk_p)(f_index_admissible, k_index)=dummy_accu_value_fk;

		}		

	}

/*(*T_fk_p).print("T_fk, line 332: after shot");*/

/*(*V_nk_p).print("V_nk, line 334: before shot");*/

	for (k_index=0; k_index<K_static; k_index++){

		for (n_index=0; n_index<N_static; n_index++){

			dummy_accu_value_kn=0;
			
			for (f_index=0; f_index<F_static_admissible; f_index++){

				f_index_admissible=f_index+1; 

				/*V_nk keeps a portion of its old amplitudes. Do not do a full update. Only a partial one, with a control scalar/parameter. Can even turn the partial update to something negligible. */
				dummy_accu_value_kn=dummy_accu_value_kn+entry_point1_fun3_uhat_KFN(k_index, f_index, n_index)/((*T_fk_p)(f_index_admissible, k_index));

				/*If you're on the last EM iteration, also update Phi_S on the way out*/
				if (em_counter == ( (( double)NUM_FEVOTTE_EM_ITERATIONS)-1 ) ) {

					/*Take a weighted average of the old phase plus the new phase, and for now depend more upon the new phase. */
					(*expj_Phi_S_nkf_p)(n_index, k_index, f_index_admissible)=0.2*(*expj_Phi_S_nkf_p)(n_index, k_index, f_index_admissible)+0.8*(entry_point1_fun3_Chat_colvecs_K_fn(k_index, f_index, n_index)/abs(entry_point1_fun3_Chat_colvecs_K_fn(k_index, f_index, n_index)));

					/*Make sure that each entry has unit amplitude. */
					(*expj_Phi_S_nkf_p)(n_index, k_index, f_index_admissible)=(*expj_Phi_S_nkf_p)(n_index, k_index, f_index_admissible)/abs((*expj_Phi_S_nkf_p)(n_index, k_index, f_index_admissible));

					(*expj_Phi_S_fkn_p)(f_index_admissible, k_index, n_index)=(*expj_Phi_S_nkf_p)(n_index, k_index, f_index_admissible);

				}

			}

			dummy_accu_value_kn=(1/((double)F_static_admissible))*dummy_accu_value_kn;

			/*(*V_nk_p)(n_index, k_index)=0.8*(*V_nk_p)(n_index, k_index)+0.2*dummy_accu_value_kn;*/
			/*(*V_nk_p)(n_index, k_index)=0.2*(*V_nk_p)(n_index, k_index)+0.8*dummy_accu_value_kn;*/
			(*V_nk_p)(n_index, k_index)=dummy_accu_value_kn;

		}

	}

/*(*V_nk_p).print("V_nk, line 370: before shot");	*/

	if  (em_counter==0){

		init_flag=false; /*Set the flag to false to signify synonymously that from this point on em_counter>0*/

	}

}	/*end EM counter*/

/*If we're out of EM and on our way back to the proposed algorithm, use some of the EM parameters to modify the Phi_S, W_fom and Z_ol parameters. */

/*entry_point1_fun3_output_phase_Phi_S_K_fn=entry_point1_fun3_Chat_colvecs_K_fn/abs(entry_point1_fun3_Chat_colvecs_K_fn);*/


}

static cube::fixed<M_static, L_static, F_static_admissible> rhat_MLF; 
static mat::fixed<M_static, L_static> rhat_ML_local; 										/*expose this*/
static cx_cube::fixed<M_static, L_static, F_static_admissible> bhat_MLF; 
static cube::fixed<M_static, L_static, F_static_admissible> bhat_MLF_real; 
static cube::fixed<M_static, L_static, F_static_admissible> bhat_MLF_imag; 

static void entry_point1_fun4_compute_rhat_bhat(void){

int f_index, l_index, m_index;

double freq_value;

/*populate rhat*/
for (l_index=0; l_index<L_static; l_index++){

	rhat_ML_local.col(l_index).zeros();

	for (f_index=0; f_index<F_static_admissible; f_index++){

		for (m_index=0; m_index<M_static; m_index++){

			freq_value=(((double)f_index+1)*((double)FS_CONSTANT))/((double)N_STFT_CONSTANT);

			rhat_MLF(m_index, l_index, f_index)=-compute_angle_tdoa_update(A_pseudoinv_EM_MLF(m_index, l_index, f_index)/A_pseudoinv_EM_MLF((int)J_VALUE, l_index, f_index))/(2*(datum::pi)*freq_value);

		}

	rhat_ML_local.col(l_index)=rhat_ML_local.col(l_index)+rhat_MLF.slice(f_index).col(l_index);

	}

	rhat_ML_local.col(l_index)=(1/((double)F_static_admissible))*rhat_ML_local.col(l_index);

}

/*populate bhat*/
for (l_index=0; l_index<L_static; l_index++){

	for (f_index=0; f_index<F_static_admissible; f_index++){

		freq_value=(((double)f_index+1)*((double)FS_CONSTANT))/((double)N_STFT_CONSTANT);

		/*bhat_MLF.slice(f_index).col(l_index)=exp(-2*(datum::pi)*freq_value*rhat_MLF.slice(f_index).col(l_index));*/

		bhat_MLF_real.slice(f_index).col(l_index)=real(cos(-2*(datum::pi)*freq_value*rhat_MLF.slice(f_index).col(l_index)));

		bhat_MLF_imag.slice(f_index).col(l_index)=imag(sin(-2*(datum::pi)*freq_value*rhat_MLF.slice(f_index).col(l_index)));

/*		bhat_MLF_real.slice(f_index).col(l_index)=real(cos(-2*(datum::pi)*freq_value*rhat_ML_local.col(l_index)));

		bhat_MLF_imag.slice(f_index).col(l_index)=imag(sin(-2*(datum::pi)*freq_value*rhat_ML_local.col(l_index)));*/

	}

}

bhat_MLF.set_real(bhat_MLF_real);
bhat_MLF.set_imag(bhat_MLF_imag);

}


static cx_cube::fixed<M_static, L_static, F_static_admissible> h_mlf_local;

static void entry_point1_fun5_compute_h_flm(mat* Z_ol_p, cube* W_fom_p, cx_cube* expj_Phi_W_fom_p){

int f_index, l_index, m_index, o_index, f_index_admissible;
cx_double dummy_value_local=0; 

for (f_index=0; f_index<F_static_admissible; f_index++){

	f_index_admissible=f_index+1;

	for (l_index=0; l_index<L_static; l_index++){

		for (m_index=0; m_index<M_static; m_index++){

			dummy_value_local=0;

			for (o_index=0; o_index<O_static; o_index++){

				dummy_value_local=dummy_value_local+((*W_fom_p)(f_index_admissible, o_index, m_index))*((*expj_Phi_W_fom_p)(f_index_admissible, o_index, m_index))*((*Z_ol_p)(o_index, l_index));

			}

			h_mlf_local(m_index, l_index, f_index)=dummy_value_local;

		}

	}

}

}

static mat::fixed<F_static, O_static> partial_d_wrt_Wfom_1;
static mat::fixed<F_static, O_static> partial_d_wrt_Wfom_2;
static cx_colvec::fixed<L_static> entry_point1_fun6_colvec1_Lx1;
static cx_colvec::fixed<L_static> entry_point1_fun6_colvec2_Lx1;

/*static mat::fixed<F_static_admissible, O_static> epsilon_mat_FO;*/

static void entry_point1_fun6_update_W_fom(mat* Z_ol_p, cube* W_fom_p, cx_cube* expj_Phi_W_fom_p){

int o_index, f_index, f_index_admissible; 

int m_index=1;

/*double gamma_local=50;*/

double gamma_local=10;

double epsilon_value; 

partial_d_wrt_Wfom_1.zeros();
partial_d_wrt_Wfom_2.zeros();

for (o_index=0; o_index<O_static; o_index++){

	for (f_index=0; f_index<F_static_admissible; f_index++){	

	f_index_admissible=f_index+1;

	entry_point1_fun6_colvec1_Lx1=conj(strans(bhat_MLF.slice(f_index).row(m_index)));

	entry_point1_fun6_colvec1_Lx1=entry_point1_fun6_colvec1_Lx1/(strans(h_mlf_local.slice(f_index).row(m_index)));

	entry_point1_fun6_colvec1_Lx1=((*expj_Phi_W_fom_p)(f_index_admissible, o_index, m_index))*entry_point1_fun6_colvec1_Lx1;

	entry_point1_fun6_colvec1_Lx1=entry_point1_fun6_colvec1_Lx1%trans(((*Z_ol_p).row(o_index)));

	partial_d_wrt_Wfom_1(f_index_admissible, o_index)=-2*as_scalar(sum(real(entry_point1_fun6_colvec1_Lx1)));

	entry_point1_fun6_colvec2_Lx1=strans(h_mlf_local.slice(f_index).row(m_index));

	entry_point1_fun6_colvec2_Lx1=entry_point1_fun6_colvec2_Lx1/strans(h_mlf_local.slice(f_index).row((int)J_VALUE));

	entry_point1_fun6_colvec2_Lx1=2*entry_point1_fun6_colvec2_Lx1/abs(entry_point1_fun6_colvec2_Lx1);

	entry_point1_fun6_colvec2_Lx1=entry_point1_fun6_colvec2_Lx1/strans(h_mlf_local.slice(f_index).row((int)J_VALUE));

	entry_point1_fun6_colvec2_Lx1=((*expj_Phi_W_fom_p)(f_index_admissible, o_index, m_index))*entry_point1_fun6_colvec2_Lx1;

	entry_point1_fun6_colvec2_Lx1=entry_point1_fun6_colvec2_Lx1%trans((*Z_ol_p).row(o_index));

	partial_d_wrt_Wfom_2(f_index_admissible, o_index)=as_scalar(sum(real(entry_point1_fun6_colvec2_Lx1)));

	epsilon_value=(1/gamma_local)*(*W_fom_p)(f_index_admissible, o_index, m_index)/partial_d_wrt_Wfom_2(f_index_admissible, o_index);

	(*W_fom_p)(f_index_admissible, o_index, m_index)=(*W_fom_p)(f_index_admissible, o_index, m_index)-epsilon_value*(partial_d_wrt_Wfom_1(f_index_admissible, o_index)+partial_d_wrt_Wfom_2(f_index_admissible, o_index));

	}

}

}

static mat::fixed<O_static, L_static> partial_d_wrt_Zol_1;
static mat::fixed<O_static, L_static> partial_d_wrt_Zol_2;
static cx_colvec::fixed<F_static_admissible> entry_point1_fun7_colvec1_Fx1;
static cx_colvec::fixed<F_static_admissible> entry_point1_fun7_colvec2_Fx1;
static cx_colvec::fixed<F_static_admissible> entry_point1_fun7_colvec3_Fx1;

static void entry_point1_fun7_update_Z_ol(mat* Z_ol_p, cube* W_fom_p, cx_cube* expj_Phi_W_fom_p){

int o_index, l_index, f_index, f_index_admissible;

int m_index=1; 

/*double gamma_local=50;*/
/*double gamma_local=10;*/
double gamma_local=1;

double epsilon_value; 

for (o_index=0; o_index<O_static; o_index++){

	for (l_index=0; l_index<L_static; l_index++){

		for (f_index=0; f_index<F_static_admissible; f_index++){

			f_index_admissible=f_index+1; 

			entry_point1_fun7_colvec1_Fx1(f_index)=conj(bhat_MLF(m_index, l_index, f_index));

			entry_point1_fun7_colvec3_Fx1(f_index)=h_mlf_local((int)J_VALUE, l_index, f_index)*((*W_fom_p)(f_index_admissible, o_index, m_index))*((*expj_Phi_W_fom_p)(f_index_admissible, o_index, m_index))-h_mlf_local(m_index, l_index, f_index)*((*W_fom_p)(f_index_admissible, o_index, (int)J_VALUE))*((*expj_Phi_W_fom_p)(f_index_admissible, o_index, (int)J_VALUE));

			entry_point1_fun7_colvec3_Fx1(f_index)=entry_point1_fun7_colvec3_Fx1(f_index)/(pow(abs(h_mlf_local(m_index, l_index, f_index)), 2));

			entry_point1_fun7_colvec1_Fx1(f_index)=entry_point1_fun7_colvec1_Fx1(f_index)*entry_point1_fun7_colvec3_Fx1(f_index);

			entry_point1_fun7_colvec2_Fx1(f_index)=h_mlf_local(m_index, l_index, f_index)/h_mlf_local((int)J_VALUE, l_index, f_index);

			entry_point1_fun7_colvec2_Fx1(f_index)=entry_point1_fun7_colvec2_Fx1(f_index)/abs(entry_point1_fun7_colvec2_Fx1(f_index));

			entry_point1_fun7_colvec2_Fx1(f_index)=entry_point1_fun7_colvec2_Fx1(f_index)*entry_point1_fun7_colvec3_Fx1(f_index);

		}

		partial_d_wrt_Zol_1(o_index, l_index)=-2*as_scalar(sum(real(entry_point1_fun7_colvec1_Fx1)));

		partial_d_wrt_Zol_2(o_index, l_index)=as_scalar(sum(real(entry_point1_fun7_colvec2_Fx1)));

		epsilon_value=(1/gamma_local)*(*Z_ol_p)(o_index, l_index)/partial_d_wrt_Zol_2(o_index, l_index);

		(*Z_ol_p)(o_index, l_index)=(*Z_ol_p)(o_index, l_index)-epsilon_value*(partial_d_wrt_Zol_1(o_index, l_index)+partial_d_wrt_Zol_2(o_index, l_index));

	}

}

}

static void entry_point1_compute_A_and_cost(cx_cube* Xtilde_fnm_p, mat* Z_ol_p, mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, cube* W_fom_p, cx_cube* expj_Phi_W_fom_p, cx_cube* expj_Phi_S_nkf_p, cx_cube* expj_Phi_S_fkn_p){

C1_Fevotte_cost_function_output_variable_colvec.print("Fevotte EM line 585:");

entry_point1_fun1_rotate_Phi_S(expj_Phi_S_nkf_p);
/*entry_point1_fun1_rotate_Phi_S(expj_Phi_S_nkf_p);*/

C1_Fevotte_cost_function_output_variable_colvec.print("Fevotte EM line 590:");

entry_point1_fun2_compute_Y_est(Y_lk_p, T_fk_p, V_nk_p);
/*entry_point1_fun2_compute_Y_est(Y_lk_p, T_fk_p, V_nk_p);*/

C1_Fevotte_cost_function_output_variable_colvec.print("Fevotte EM line 595:");

entry_point1_fun0_rotate_X_Y(Xtilde_fnm_p);
/*entry_point1_fun3_fun0_rotate_X_Y(Xtilde_fnm_p);*/

entry_point1_fun3_compute_A_via_EM(Y_lk_p, T_fk_p, V_nk_p, expj_Phi_S_nkf_p, expj_Phi_S_fkn_p);
/*entry_point1_fun3_fun1_compute_A_via_EM();*/

(*(expj_Phi_S_nkf_p)).elem( find_nonfinite((*(expj_Phi_S_nkf_p))) ).fill(1);
(*(expj_Phi_S_fkn_p)).elem( find_nonfinite((*(expj_Phi_S_fkn_p))) ).fill(1);

/*Print T and V for debug*/
/*(*T_fk_p).print("T_fk, line 613 whatup:");
(*V_nk_p).print("V_nk, line 614 whatup:");
(*expj_Phi_S_nkf_p).slice(0).print("expj_Phi_S_nkf_p, line 615 whatup home slice(0)");
(*expj_Phi_S_nkf_p).slice(1).print("expj_Phi_S_nkf_p, line 615 whatup home slice(1)");*/

C1_Fevotte_cost_function_output_variable_colvec.print("Fevotte EM line 600:");

/*entry_point1_fun4_compute_rhat_bhat();*/
/*entry_point1_fun3_fun1_compute_rhat_bhat();*/

C1_Fevotte_cost_function_output_variable_colvec.print("Fevotte EM line 605:");

/*entry_point1_fun5_compute_h_flm(Z_ol_p, W_fom_p, expj_Phi_W_fom_p);*/
/*entry_point1_fun4_compute_h_flm(Z_ol_p, W_fom_p, expj_Phi_W_fom_p);*/

C1_Fevotte_cost_function_output_variable_colvec.print("Fevotte EM line 610:");

/*entry_point1_fun6_update_W_fom(Z_ol_p, W_fom_p, expj_Phi_W_fom_p);*/

/*entry_point1_fun7_update_Z_ol(Z_ol_p, W_fom_p, expj_Phi_W_fom_p);*/
/*entry_point1_fun6_update_Z_ol(Z_ol_p, W_fom_p, expj_Phi_W_fom_p);*/

C1_Fevotte_cost_function_output_variable_colvec.print("Fevotte EM line 615:");

}

void Fevotte_EM_update_module1_entry_point1(cx_cube* Xtilde_fnm_p, mat* Z_ol_p, mat* Y_lk_p, mat* T_fk_p, mat* V_nk_p, cube* W_fom_p, cx_cube* expj_Phi_W_fom_p, cx_cube* expj_Phi_S_nkf_p, cx_cube* expj_Phi_S_fkn_p){

entry_point1_compute_A_and_cost(Xtilde_fnm_p, Z_ol_p, Y_lk_p, T_fk_p, V_nk_p, W_fom_p, expj_Phi_W_fom_p, expj_Phi_S_nkf_p, expj_Phi_S_fkn_p);
	
}